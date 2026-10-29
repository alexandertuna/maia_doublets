#
# Set up the MuColl software stack and the maia_doublets venv.
#
# Use interactively:            source run/env.sh
# or from scripts / apptainer:  source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
#
# The container's setup_mucoll is an alias from /etc/profile.d/aliases.sh,
# which only exists in interactive login shells. So source its target directly.
#
RUN_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source /opt/spack/opt/spack/__spack_path_placeholder__/__spack_path_placeholder__/__spack_path_placeholder__/__spack_path_placeholder__/linux-x86_64/mucoll-stack-*/setup.sh
source ${RUN_DIR}/env_ubuntu24/bin/activate
