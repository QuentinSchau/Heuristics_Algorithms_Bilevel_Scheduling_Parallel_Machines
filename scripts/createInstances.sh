#!/bin/bash

# Define the range of values for the variables
seed=0
nbInstGen=10
n_values=(40 50 60 70 80)
frac_n=(0.25 0.5 0.75)
idJobs=0

# Get the absolute directory of the script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Generation of instances
# Loop over each value of n
for n in "${n_values[@]}"; do
    python3 $(realpath $SCRIPT_DIR/generateInstances.py) \
              --no-identical-jobs $idJobs \
              --frac-of-n "${frac_n[@]}" \
              --m 1,1 --seed $seed --nb-instance-to-generate $nbInstGen \
              --n-list $n --config-file-name-generate config_generate_N${n}_M2 \
              --path-save-instance /instances/N$n/instances/

    $(realpath $SCRIPT_DIR/../bin/Release/bilevel-scheduling) \
      $(realpath $SCRIPT_DIR/../instances/configs/config_generate_N${n}_M2.json)

    # Second instance generation and execution
    python3 $(realpath $SCRIPT_DIR/generateInstances.py) \
      --no-identical-jobs $idJobs \
      --frac-of-n "${frac_n[@]}" \
      --m 2,2 --seed $seed --nb-instance-to-generate $nbInstGen \
      --n-list $n --config-file-name-generate config_generate_N${n}_M4 \
      --path-save-instance /instances/N$n/instances/

    $(realpath $SCRIPT_DIR/../bin/Release/bilevel-scheduling) \
      $(realpath $SCRIPT_DIR/../instances/configs/config_generate_N${n}_M4.json)

    echo "Completed generation of instances"
done
