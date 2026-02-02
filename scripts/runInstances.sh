#!/bin/bash

# Number of parallel jobs (set to the number of CPUs or desired concurrency)

# Get the absolute directory of the script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

num_cpus=6
# Default behavior: Use the provided find method
find "$(realpath "$SCRIPT_DIR/../instances/configs/config_solve/")" -type f -name "*.json" | \
parallel -j "$num_cpus" --bar --joblog execution.log --resume \
"$(realpath "$SCRIPT_DIR/../bin/Release/bilevel-scheduling") {} >> bilevel_scheduling.log 2>&1"


