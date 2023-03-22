#!/bin/bash

export OMP_NUM_THREADS=1
export cur_dir=`pwd`
export main_dir="/home/kkiesling/depletion/hybrid-depletion"

function run_direct () {
    mkdir -p ${cur_dir}/direct/
    cd ${cur_dir}/direct/
    ln -s ${main_dir}/model/chain_endfb71_pwr.xml .
    python ${main_dir}/deplete/reduce_chain.py -c ${main_dir}/model/chain_endfb71_pwr.xml
    echo "Running Direct"
    mpiexec -n 20 python ${main_dir}/deplete/run_depletion.py -d -m ${main_dir}/model/pin/ -c ${main_dir}/model/chain_endfb71_pwr_reduced.xml
}

function run_flux () {
    for g in 300 500 2500 10000
    do
        mkdir -p ${cur_dir}/flux/${g}
        cd ${cur_dir}/flux/${g}
        ln -s ${main_dir}/model/chain_endfb71_pwr.xml .
        python ${main_dir}/deplete/reduce_chain.py -c ${main_dir}/model/chain_endfb71_pwr.xml
        echo "Running Flux ${g}"
        mpiexec -n 20 python ${main_dir}/deplete/run_depletion.py -f -g ${g} -m ${main_dir}/model/pin/ -c ${main_dir}/model/chain_endfb71_pwr_reduced.xml
    done
}

run_direct
run_flux