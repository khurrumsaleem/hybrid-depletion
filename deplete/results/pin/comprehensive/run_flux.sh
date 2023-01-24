#!/bin/bash

export OMP_NUM_THREADS=1
export cur_dir=`pwd`
export main_dir="/home/kkiesling/depletion/hybrid-depletion"

function run_calcs () {
    groups=$@
    for g in ${groups[@]}
        do
            cd ${cur_dir}/flux/${g}
            cp ${main_dir}/model/chain_endfb71_pwr.xml .
            echo "Running Flux ${g}"
            mpiexec -n 20 python ${main_dir}/deplete/run_depletion.py -g ${g} -f -m ${main_dir}/model/pin/ -c chain_endfb71_pwr.xml
        done
}

g=(10000)
run_calcs ${g[@]}
