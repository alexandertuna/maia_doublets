#!/usr/bin/env bash
#
# Scan the C++ maia_doublets timing on v06 or v07 BIB samples (MAIA_DETECTOR, default v06),
# running cpp/run.sh once per
# (time window, file). Each BIB file has one event.
#
# Grid: every BIB fraction (10-120%) x time window {maia, wide, none}.
# Both sensors of each double layer are always kept (the algorithm needs them).
# 110% and 120% are fine here even without vertex-barrel hits: only IT/OT hits are used.
#
# Run on the host, not inside the container: each job starts its own container,
# so that /usr/bin/time can record its peak memory. Made for a screen session:
# > screen -S doublets_scan
# > ./cpp/run_scan.sh
# > MAIA_SCAN_FILES=3 MAIA_SCAN_JOBS=4 ./cpp/run_scan.sh
# > MAIA_DETECTOR=v07 ./cpp/run_scan.sh
#
# Jobs whose log already has the timing table are skipped, so the scan can be
# stopped (ctrl-c) and restarted, or extended by raising MAIA_SCAN_FILES.
# Per job: output/${GEO}_${DATASET}_digi_${SMEAR}/doublets_*.log (timing) and doublets_*.time (peak memory).
# Progress: output/doublets_scan.log
#
set -eo pipefail
cd "$(dirname "$0")"

FILES=${MAIA_SCAN_FILES:-1}
JOBS=${MAIA_SCAN_JOBS:-4}
GEO=${MAIA_DETECTOR:-v06}
SMEAR=${MAIA_RESOLUTION:-10um}

CONTAINER=/cvmfs/unpacked.cern.ch/ghcr.io/muoncollidersoft/mucoll-sim-ubuntu24:v2.9.8-amd64
SCAN_LOG=../output/doublets_scan.log

# Biggest first.
DATASETS="neutrinoGun120 neutrinoGun110 neutrinoGun100 neutrinoGun90 neutrinoGun80 neutrinoGun70
          neutrinoGun60 neutrinoGun50 neutrinoGun40 neutrinoGun30 neutrinoGun20 neutrinoGun10"
WINDOWS="none wide maia"

# One line per job: DATASET WINDOW FILE
list_jobs() {
  for FILE in $(seq 0 $((FILES - 1))); do
    for DATASET in ${DATASETS}; do
      for WINDOW in ${WINDOWS}; do
        echo ${DATASET} ${WINDOW} ${FILE}
      done
    done
  done
}

run_job() {
  read DATASET WINDOW FILE <<< "$1"
  # Same naming as run.sh
  NAME=../output/${GEO}_${DATASET}_digi_${SMEAR}/doublets_${WINDOW}_${FILE}
  if grep -q "AlgTimingAuditor .*EVENT LOOP" ${NAME}.log 2>/dev/null; then
    echo "$(date '+%F %T') skip  ${GEO} $1" | tee -a ${SCAN_LOG}
    return 0
  fi
  echo "$(date '+%F %T') start ${GEO} $1" | tee -a ${SCAN_LOG}
  mkdir -p $(dirname ${NAME})
  # 128 GB virtual-memory cap per job, to protect the shared node
  if ( ulimit -v 134217728
       /usr/bin/time -v -o ${NAME}.time \
         apptainer exec ${CONTAINER} env \
           MAIA_DETECTOR=${GEO} \
           MAIA_RESOLUTION=${SMEAR} \
           MAIA_DATASET=${DATASET} \
           MAIA_FILE=${FILE} \
           MAIA_TIME_WINDOW=${WINDOW} \
           ./run.sh > /dev/null 2>&1 ); then
    echo "$(date '+%F %T') done  ${GEO} $1" | tee -a ${SCAN_LOG}
  else
    echo "$(date '+%F %T') FAIL  ${GEO} $1" | tee -a ${SCAN_LOG}
  fi
}
export -f run_job
export GEO SMEAR CONTAINER SCAN_LOG

mkdir -p ../output
echo "$(date '+%F %T') scan ${GEO}: ${FILES} file(s) per point, ${JOBS} parallel jobs, $(list_jobs | wc -l) jobs" | tee -a ${SCAN_LOG}
list_jobs | xargs -P ${JOBS} -I{} bash -c 'run_job "{}"'
echo "$(date '+%F %T') scan finished" | tee -a ${SCAN_LOG}
