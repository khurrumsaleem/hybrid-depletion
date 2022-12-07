#!/bin/bash

export OMP_NUM_THREADS=20
export cur_dir=`pwd`
export main_dir="/home/kkiesling/depletion/hybrid-depletion"

function run_calcs () {
    for h in 1 2
    do
        for n in all actinides mix
        do
            for g in 300 500 2500 10000
            do
                cd ${cur_dir}/hybrid${h}/${n}/${g}
                echo "Running Hybrid ${h} ${n} ${g}"
                mpiexec -n 2 python ${main_dir}/deplete/run_depletion.py -m ${main_dir}/model/pin/ -c ${main_dir}/model/chain_endfb71_pwr.xml -n ${n} -g ${g} -y ${h}
            done
        done
    done
}

run_calcs