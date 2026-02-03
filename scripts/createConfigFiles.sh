#!/bin/bash

#
# Copyright (C) 2024
# Laboratoire d'Informatique Fondamentale et Appliquée de Tours, Tours, France
#
# DIGEP, Politecnico di Torino, Corso Duca degli Abruzzi 24, Torino, Italy
# This file is part of bilevel-scheduling.
#
# bilevel-scheduling is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.
#
# bilevel-scheduling is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with bilevel-scheduling. If not, see <https://www.gnu.org/licenses/>.
#

# Define the range of values for the variables
nbInstGen=10
methods=("BeamSearch")
strategies=("local-search")
alphas=(0.0 0.1)
versions=(1 3)
versionLS=1
bestImprove=0
n_values=(40 50 60 70 80)
frac_n=(0.25 0.5 0.75)
timeLimit=60
timeLimitLS=300
beamSizes=(1)
recovering=(1)
nbSolutionForMSLS=1000
autoSetting=1 # use auto setting in MSLS
# Get the absolute directory of the script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Generation of instances
# Loop over each method to generate json file configuration
for n in "${n_values[@]}"; do
  for method in "${methods[@]}"; do
    # Check if method is "BeamSearch"
    if [[ "$method" == "BeamSearch" ]]; then
      for reco in "${recovering[@]}"; do       
        for beamSize in "${beamSizes[@]}"; do
          for alpha in "${alphas[@]}"; do
            if [[ "$reco" == "1" ]]; then
              # Loop over each strategy
              for strategy in "${strategies[@]}"; do
                if [[ "$strategy" == "local-search" ]]; then
                  for version in "${versions[@]}"; do
                    echo "Running with method=$method, n=$n, strategy=$strategy, alpha= $alpha,nbSolutionForMSLS=$nbSolutionForMSLS, use reco=$reco, version= $version and beam size=$beamSize"

                    # First instance generation and execution
                    python3 "$(realpath $SCRIPT_DIR/generateInstances.py)" \
                      --frac-of-n "${frac_n[@]}" --timeLimit $timeLimit --timeLimitLS $timeLimitLS \
                      --m 1,1 --nb-instance-to-generate $nbInstGen \
                      --n-list $n --output-result /instances/N$n/2M_ --methods $method \
                      --config-file-name-solve config_solve/config_solve_N${n}_M2_${method}_beamSize_${beamSize}_alpha_${alpha}_reco_${reco}_nbSolutionForMSLS_${nbSolutionForMSLS}_autoSetting_${autoSetting}_recoStrat_${strategy}_versionRBS_${version}_versionLS_${versionLS}_bestImp_${bestImprove}_ \
                      --beamSize $beamSize --alpha $alpha --versionRBS $version --versionLS $versionLS --bestImproveLS $bestImprove --recovering $reco --strategy $strategy --alpha $alpha --nbSolutionForMSLS $nbSolutionForMSLS --autoSetting $autoSetting \
                      --path-save-instance /instances/N$n/instances/

                    # Second instance generation and execution
                    python3 "$(realpath $SCRIPT_DIR/generateInstances.py)" \
                      --frac-of-n "${frac_n[@]}" --timeLimit $timeLimit --timeLimitLS $timeLimitLS \
                      --m 2,2 --nb-instance-to-generate $nbInstGen \
                      --n-list $n --output-result /instances/N$n/4M_ --methods $method \
                      --config-file-name-solve config_solve/config_solve_N${n}_M4_${method}_beamSize_${beamSize}_alpha_${alpha}_reco_${reco}_nbSolutionForMSLS_${nbSolutionForMSLS}_autoSetting_${autoSetting}_recoStrat_${strategy}_versionRBS_${version}_versionLS_${versionLS}_bestImp_${bestImprove}_ \
                      --beamSize $beamSize --alpha $alpha --versionRBS $version --versionLS $versionLS --bestImproveLS $bestImprove --recovering $reco --strategy $strategy --alpha $alpha --nbSolutionForMSLS $nbSolutionForMSLS --autoSetting $autoSetting \
                      --path-save-instance /instances/N$n/instances/
                  done
                else
                  echo "Running with method=$method, n=$n, strategy=$strategy, alpha= $alpha,nbSolutionForMSLS=$nbSolutionForMSLS, use reco=$reco and beam size=$beamSize"

                  # First instance generation and execution
                  python3 "$(realpath $SCRIPT_DIR/generateInstances.py)" \
                    --frac-of-n "${frac_n[@]}" --timeLimit $timeLimit --timeLimitLS $timeLimitLS \
                    --m 1,1 --nb-instance-to-generate $nbInstGen \
                    --n-list $n --output-result /instances/N$n/2M_ --methods $method \
                    --config-file-name-solve config_solve/config_solve_N${n}_M2_${method}_beamSize_${beamSize}_alpha_${alpha}_reco_${reco}_nbSolutionForMSLS_${nbSolutionForMSLS}_autoSetting_${autoSetting}_recoStrat_${strategy}_ \
                    --beamSize $beamSize --alpha $alpha --recovering $reco --strategy $strategy --nbSolutionForMSLS $nbSolutionForMSLS --autoSetting $autoSetting \
                    --path-save-instance /instances/N$n/instances/

                  # Second instance generation and execution
                  python3 "$(realpath $SCRIPT_DIR/generateInstances.py)" \
                    --frac-of-n "${frac_n[@]}" --timeLimit $timeLimit --timeLimitLS $timeLimitLS \
                    --m 2,2 --nb-instance-to-generate $nbInstGen \
                    --n-list $n --output-result /instances/N$n/4M_ --methods $method \
                    --config-file-name-solve config_solve/config_solve_N${n}_M4_${method}_beamSize_${beamSize}_alpha_${alpha}_reco_${reco}_nbSolutionForMSLS_${nbSolutionForMSLS}_autoSetting_${autoSetting}_recoStrat_${strategy}_ \
                    --beamSize $beamSize --alpha $alpha --recovering $reco --strategy $strategy --nbSolutionForMSLS $nbSolutionForMSLS --autoSetting $autoSetting \
                    --path-save-instance /instances/N$n/instances/
                fi
                echo "Completed for method=$method, n=$n, reco=$reco,strategy=$strategy, alpha= $alpha,nbSolutionForMSLS=$nbSolutionForMSLS, and beam size=$beamSize"
              done
            else
              echo "Running with method=$method, n=$n, use reco=$reco, alpha= $alpha,nbSolutionForMSLS=$nbSolutionForMSLS, and beam size=$beamSize"

              # First instance generation and execution
              python3 "$(realpath $SCRIPT_DIR/generateInstances.py)" \
                --frac-of-n "${frac_n[@]}" --timeLimit $timeLimit --timeLimitLS $timeLimitLS \
                --m 1,1 --nb-instance-to-generate $nbInstGen \
                --n-list $n --output-result /instances/N$n/2M_ --methods $method \
                --config-file-name-solve config_solve/config_solve_N${n}_M2_${method}_beamSize_${beamSize}_alpha_${alpha}_reco_${reco}_nbSolutionForMSLS_${nbSolutionForMSLS}_autoSetting_${autoSetting}_versionLS_${versionLS}_bestImp_${bestImprove}_ \
                --beamSize $beamSize --alpha $alpha --recovering $reco --versionLS $versionLS --bestImproveLS $bestImprove --nbSolutionForMSLS $nbSolutionForMSLS --autoSetting $autoSetting \
                --path-save-instance /instances/N$n/instances/

              # Second instance generation and execution
              python3 "$(realpath $SCRIPT_DIR/generateInstances.py)" \
                --frac-of-n "${frac_n[@]}" --timeLimit $timeLimit --timeLimitLS $timeLimitLS \
                --m 2,2 --nb-instance-to-generate $nbInstGen \
                --n-list $n --output-result /instances/N$n/4M_ --methods $method \
                --config-file-name-solve config_solve/config_solve_N${n}_M4_${method}_beamSize_${beamSize}_alpha_${alpha}_reco_${reco}_nbSolutionForMSLS_${nbSolutionForMSLS}_autoSetting_${autoSetting}_versionLS_${versionLS}_bestImp_${bestImprove}_ \
                --beamSize $beamSize --alpha $alpha --recovering $reco --versionLS $versionLS --bestImproveLS $bestImprove --nbSolutionForMSLS $nbSolutionForMSLS --autoSetting $autoSetting \
                --path-save-instance /instances/N$n/instances/
              echo "Completed for method=$method, n=$n, reco=$reco,nbSolutionForMSLS=$nbSolutionForMSLS, and beam size=$beamSize"
            fi
          done
        done
      done
    else
      for version in "${versions[@]}"; do
        # If method is other heuristic
        echo "Running with method=$method, n=$n and version=$version"

        # First instance generation and execution
        python3 "$(realpath $SCRIPT_DIR/generateInstances.py)" \
          --frac-of-n "${frac_n[@]}" --timeLimit $timeLimit \
          --m 1,1 --nb-instance-to-generate $nbInstGen \
          --n-list $n --versionLS $version --bestImproveLS $bestImprove --output-result /instances/N$n/2M_ --methods $method \
          --config-file-name-solve config_solve/config_solve_N${n}_M2_${method}_versionLS_${version}_bestImp_${bestImprove}_ \
          --path-save-instance /instances/N$n/instances/

        # Second instance generation and execution
        python3 "$(realpath $SCRIPT_DIR/generateInstances.py)" \
          --frac-of-n "${frac_n[@]}" --timeLimit $timeLimit \
          --m 2,2 --nb-instance-to-generate $nbInstGen \
          --n-list $n --versionLS $version --bestImproveLS $bestImprove --output-result /instances/N$n/4M_ --methods $method \
          --config-file-name-solve config_solve/config_solve_N${n}_M4_${method}_versionLS_${version}_bestImp_${bestImprove}_ \
          --path-save-instance /instances/N$n/instances/
        echo "Completed for method=$method and n=$n"
      done
    fi
  done
done
