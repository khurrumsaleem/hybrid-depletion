export OMP_NUM_THREADS=1
export cur_dir=`pwd`
export main_dir="/home/kkiesling/depletion/hybrid-depletion"

function run_hybrid () {
    for g in 300 500 2500 10000
    do
        for n in actinides mix
        do
            mkdir -p ${cur_dir}/hybrid2/${n}/${g}
            cd ${cur_dir}/hybrid2/${n}/${g}
            ln -s ${main_dir}/model/chain_endfb71_pwr.xml .
            python ${main_dir}/deplete/reduce_chain.py -c ${main_dir}/model/chain_endfb71_pwr.xml
            echo "Running Hybrid 2 ${g} ${n}"
            mpiexec -n 20 python ${main_dir}/deplete/run_depletion.py -y 2 -n ${n} -g ${g} -m ${main_dir}/model/pin/ -c ${main_dir}/model/chain_endfb71_pwr_reduced.xml
        done
    done
}

run_hybrid
