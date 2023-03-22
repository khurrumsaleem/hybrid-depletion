import openmc
import openmc.deplete
import argparse
import pathlib

def parse_args():
    ###############################################################################
    #                           Parse Args
    ###############################################################################
    parser = argparse.ArgumentParser(description='reduce chain file')
    parser.add_argument('-c', '--chain', type=pathlib.Path,
                        help='Path to chain .xml file and chain file to load.')
    args = parser.parse_args()

    return args


if __name__ == "__main__":
    args = parse_args()
    chain_file = str(args.chain)
    chain_red_path = chain_file[:-4] + '_reduced.xml'

    chain_full = openmc.deplete.Chain.from_xml(chain_file)
    stable = [
        nuc.name
        for nuc in chain_full.nuclides
        if nuc.half_life is None or nuc.half_life > 1e15
    ]
    chain_reduced = chain_full.reduce(stable)
    chain_reduced.export_to_xml(chain_red_path)
