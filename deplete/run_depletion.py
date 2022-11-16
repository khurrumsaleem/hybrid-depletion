import openmc
import openmc.deplete
import numpy as np
import shutil
import argparse
import pathlib
import os

openmc.deplete.pool.NUM_PROCESSES = 20
openmc.deplete.pool.USE_MULTIPROCESSING = True

def parse_args():
    ###############################################################################
    #                           Parse Args
    ###############################################################################
    parser = argparse.ArgumentParser(description='Specify Simulations to run')
    parser.add_argument('-g', '--groups', default=500, type=int, choices=[300, 500],
                        help='Number of energy groups to use; 500 or 300 (default=500)')
    parser.add_argument('-d', '--direct', action='store_true',
                        help='Run direct tally')
    parser.add_argument('-f', '--flux', action='store_true',
                        help='Run full flux tally')
    parser.add_argument('-y', '--hybrid', type=int, choices=[1, 2],
                        help='Hybrid tally with either 1 or 2 direct reaction rates')
    parser.add_argument('-a', '--all', action='store_true',
                        help='Run all types of depletion calcs')
    parser.add_argument('-m', '--model', default='.', type=pathlib.Path,
                        help='Path to model .xml files and chain file to load. Default cwd.')
    parser.add_argument('-c', '--chain', type=pathlib.Path,
                        help='Path to chain .xml file and chain file to load.')
    args = parser.parse_args()

    return args


if __name__ == "__main__":
    args = parse_args()

    ###############################################################################
    #                   Load Model (../make_pin_model.py)
    ###############################################################################
    model_path = str(args.model)
    model = openmc.model.Model.from_xml(geometry=model_path + '/geometry.xml',
                                        settings=model_path + '/settings.xml',
                                        materials=model_path + '/materials.xml')

    ###############################################################################
    #                 Energy Group Structure (Salcedo-Perez 2019 M&C)
    ###############################################################################
    groups500 = list(np.logspace(np.log10(1e-5), np.log10(400e3), 50, endpoint=False)) +\
                list(np.logspace(np.log10(400e3), np.log10(3.21e6), 100, endpoint=False)) +\
                list(np.logspace(np.log10(3.21e6), np.log10(8.025e6), 260)) +\
                list(np.logspace(np.log10(8.05e6), np.log10(20e6), 90))
    groups300 = list(np.logspace(np.log10(1e-5), np.log10(400e3), 10, endpoint=False)) +\
                list(np.logspace(np.log10(400e3), np.log10(3.21e6), 10, endpoint=False)) +\
                list(np.logspace(np.log10(3.21e6), np.log10(8.025e6), 260)) +\
                list(np.logspace(np.log10(8.05e6), np.log10(20e6), 20))

    if args.groups == 300:
        groups = groups300
    else:
        # default / other option
        groups = groups500

    ###############################################################################
    #                  Reduce Chain File
    ###############################################################################
    chain_file = str(args.chain)
    chain_red_path = chain_file[:-4] + '_reduced.xml'
    #chain_file = chain_path + '/chain_endfb71_pwr.xml'

    # don't reduce if already reduced
    if not os.path.exists(chain_red_path):
        chain_full = openmc.deplete.Chain.from_xml(chain_file)
        stable = [
            nuc.name
            for nuc in chain_full.nuclides
            if nuc.half_life is None or nuc.half_life > 1e15
        ]
        chain_reduced = chain_full.reduce(stable)
        chain_reduced.export_to_xml(chain_red_path)
        chain_file = chain_red_path

    ###############################################################################
    #                  Set burnup steps
    ###############################################################################
    # cumulative steps in MWd/kg
    burnup_cum = np.array([
        0.1, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0,
        12.5, 15.0, 17.5, 20.0, 22.5, 25.0, 27.5, 30.0, 32.5, 35.0, 37.5,
        40.0, 42.5, 45.0, 47.5, 50.0
    ])
    burnup = np.diff(burnup_cum, prepend=0.0)
    power = 174  # W/cm

    ###############################################################################
    #             Run depletion with all nuclides directly tallied
    ###############################################################################
    if args.all or args.direct:
        print("\n******* Running Direct Calculation *******\n")
        op = openmc.deplete.CoupledOperator(model, chain_file=chain_file,
                                            reaction_rate_mode='direct')
        integrator = openmc.deplete.PredictorIntegrator(op, burnup, power,
                                                        timestep_units='MWd/kg')
        integrator.integrate()

    ###############################################################################
    #                Run depletion with all nuclides flux tallied
    ###############################################################################
    if args.all or args.flux:
        print("\n******* Running Flux Calculation *******\n")
        op = openmc.deplete.CoupledOperator(model, chain_file=chain_file,
                                            reaction_rate_mode='flux',
                                            reaction_rate_opts={'energies': groups,
                                            'reactions': None})
        integrator = openmc.deplete.PredictorIntegrator(op, burnup, power,
                                                        timestep_units='MWd/kg')
        integrator.integrate()

    ###############################################################################
    #                          Run hybrid depletion - both reaction rates
    ###############################################################################
    rr2 = ['fission', '(n,gamma)']
    if args.all or args.hybrid == 2:
        print("\n******* Running Hybrid 2 Calculation *******\n")
        # should direct tally all nuclides for each reaction rate
        op = openmc.deplete.CoupledOperator(model, chain_file=chain_file,
                                            reaction_rate_mode='flux',
                                            reaction_rate_opts={'energies': groups,
                                            'reactions': rr2})
        integrator = openmc.deplete.PredictorIntegrator(op, burnup, power,
                                                        timestep_units='MWd/kg')
        integrator.integrate()

    ###############################################################################
    #                          Run hybrid depletion - one reaction rates
    ###############################################################################
    rr1 = ['(n,gamma)']
    if args.all or args.hybrid == 1:
        print("\n******* Running Hybrid 1 Calculation *******\n")
        op = openmc.deplete.CoupledOperator(model, chain_file=chain_file,
                                            reaction_rate_mode='flux',
                                            reaction_rate_opts={'energies': groups,
                                            'reactions': rr1})
        integrator = openmc.deplete.PredictorIntegrator(op, burnup, power,
                                                        timestep_units='MWd/kg')
        integrator.integrate()
