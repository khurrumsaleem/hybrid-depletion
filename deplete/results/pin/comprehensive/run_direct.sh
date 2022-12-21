#!/bin/bash

export OMP_NUM_THREADS=1
export cur_dir=`pwd`
export main_dir="/home/kkiesling/depletion/hybrid-depletion"

function run_calcs () {
    for g in 300 500 2500 10000
    do
        cd ${cur_dir}/direct/${g}
        echo "Running Direct ${g}"
        mpiexec -n 20 python ${main_dir}/deplete/run_depletion.py -m ${main_dir}/model/pin/ -c ${main_dir}/model/chain_endfb71_pwr.xml -g ${g} -d
    done
}

run_calcs