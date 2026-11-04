#
# Add the C++ port's Marlin processor to MARLIN_DLL. Source after run/env.sh:
# > source cpp/env.sh
#
export MAIA_DOUBLETS_CPP_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export MARLIN_DLL=${MARLIN_DLL}:${MAIA_DOUBLETS_CPP_DIR}/install/lib/libMaiaDoublets.so
