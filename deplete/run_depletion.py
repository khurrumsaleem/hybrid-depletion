import openmc
import openmc.deplete
import numpy as np
import argparse
import pathlib
import os

#openmc.deplete.pool.NUM_PROCESSES = 20
#openmc.deplete.pool.USE_MULTIPROCESSING = True

def parse_args():
    ###############################################################################
    #                           Parse Args
    ###############################################################################
    parser = argparse.ArgumentParser(description='Specify Simulations to run')
    parser.add_argument('-g', '--groups', default=500, type=int, choices=[300, 500, 2500, 10000],
                        help='Number of energy groups to use; 500 or 300 or 2500 (default=500)')
    parser.add_argument('-n', '--nuclides', default='all', type=str, choices=['all', 'actinides', 'mix'],
                        help="Which nuclides to direct tally")
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
    # Energy Group Structure (Salcedo-Perez 2019 M&C)/ECP Milestone report 10/1/2020
    ###############################################################################
    if args.groups == 10000:
        groups = list(np.logspace(np.log10(1e-5), np.log10(2e7), 10000 + 1))
    else:
        if args.groups == 300:
            ng = [10, 10, 260, 20]
        elif args.groups == 500:
            ng = [10, 100, 260, 90]
        elif args.groups == 2500:
            ng = [250, 500, 1300, 450]

        groups = list(np.logspace(np.log10(1e-5), np.log10(4e5), ng[0], endpoint=False)) +\
                list(np.logspace(np.log10(4e5), np.log10(3.21e6), ng[1], endpoint=False)) +\
                list(np.logspace(np.log10(3.21e6), np.log10(8.025e6), ng[2], endpoint=False)) +\
                list(np.logspace(np.log10(8.025e6), np.log10(2e7), ng[3] + 1))

    ###############################################################################
    #                  Nuclides to direct tally
    ##############################################################################
    if args.nuclides == 'all':
        nuclides = None
    elif args.nuclides == 'actinides':
        nuclides = ['Th230', 'Th231', 'Th232', 'Th234', 'Pa231', 'Pa232',
                    'Pa233', 'U232', 'U233', 'U234', 'U235', 'U236', 'U237',
                    'U238', 'U239', 'Np236', 'Np237', 'Np238', 'Np239',
                    'Pu236', 'Pu238', 'Pu239', 'Pu240', 'Pu241', 'Pu242',
                    'Pu243', 'Am241', 'Am242', 'Am242_m1', 'Am243', 'Cm242', 'Cm244', 'Xe135']
    elif args.nuclides == 'mix':
        nuclides = ['Kr83', 'Tc99', 'Ru101', 'Rh103', 'Rh105', 'Xe131', 'Xe133',
                    'Xe135', 'Cs133', 'Cs135', 'Pr143', 'Nd143', 'Nd145',
                    'Pm147', 'Pm148_m1', 'Pm149', 'Pm151', 'Sm149', 'Sm150',
                    'Sm151', 'Sm152', 'Sm153', 'Eu153', 'Eu155', 'Gd157', 'U234',
                    'U235', 'U236', 'U237', 'U238', 'Pu239', 'Pu240', 'U239', 'Np239']

    print("Tallying these direct nuclides:")
    print(nuclides)

    ###############################################################################
    #                  Reduce Chain File
    ###############################################################################
    chain_file = str(args.chain)

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
                                            'reactions': rr2, 'nuclides': nuclides})
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
                                            'reactions': rr1, 'nuclides': nuclides})
        integrator = openmc.deplete.PredictorIntegrator(op, burnup, power,
                                                        timestep_units='MWd/kg')
        integrator.integrate()
