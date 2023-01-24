#!/bin/bash

export OMP_NUM_THREADS=1
export cur_dir=`pwd`
export main_dir="/home/kkiesling/depletion/hybrid-depletion"

function run_calcs () {
    n=$1
    shift
    groups=$@
    h=2
    for g in ${groups[@]}
    do
        cd ${cur_dir}/hybrid${h}/${n}/${g}
        cp ${main_dir}/model/chain_endfb71_pwr.xml .
        echo "Running Hybrid 2 ${n} ${g}"
        mpiexec -n 20 python ${main_dir}/deplete/run_depletion.py -n ${n} -g ${g} -y ${h} -m ${main_dir}/model/pin/ -c chain_endfb71_pwr.xml
    done
}

gr=(300 500 2500 10000)
run_calcs all ${gr[@]}
run_calcs actinides ${gr[@]}
run_calcs mix ${gr[@]}
