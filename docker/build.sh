#!/bin/bash
set -euo pipefail

for docker_path in $(find . -name Dockerfile); do pushd $(dirname $docker_path); docker build -t $(basename ${PWD}):latest .; popd; done


