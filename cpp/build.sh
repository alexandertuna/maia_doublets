#!/usr/bin/env bash
#
# Build the C++ port into cpp/build and install it into cpp/install.
#
# Run inside the container:
# > ./cpp/build.sh
#
set -eo pipefail
cd "$(dirname "$0")"
source ../run/env.sh

set -x
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=install
cmake --build build -j8
cmake --install build
