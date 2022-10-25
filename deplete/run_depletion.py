import openmc
import openmc.deplete
import numpy as np
import shutil

openmc.deplete.pool.NUM_PROCESSES = 1
openmc.deplete.pool.USE_MULTIPROCESSING = False

def setup_flux_operator(reactions, nuclides):
    op = openmc.deplete.CoupledOperator(model, chain_file=chain_file,
                                        fission_yield_mode="average",
                                        reaction_rate_mode='flux',
                                        reaction_rate_opts={'energies': groups500,
                                                            'reactions': reactions,
                                                            'nuclides': nuclides})
    return op

def run(op):
    integrator = openmc.deplete.PredictorIntegrator(op, burnup, power,
                                                    timestep_units='MWd/kg')
    integrator.integrate()

def mv_results(dest):
    shutil.move("depletion_results.h5", dest)

###############################################################################
#                   Load Model (../make_pin_model.py)
###############################################################################
model = openmc.model.Model.from_xml()

###############################################################################
#                   List of Nuclides (Romano 2021)
###############################################################################
actinides = ['U234', 'U235', 'U236', 'U238', 'U239','Np239',
             'Pu238', 'Pu239', 'Pu240', 'Pu241', 'Pu242',
             'Am241', 'Am242', 'Am242_m1', 'Am243', 'Am244',
             'Cm242', 'Cm243', 'Cm244', 'Cm245', 'Cm246']
fps = ['Kr85', 'Sr90', 'Y90', 'Zr93', 'Mo95', 'Mo97', 'Tc99', 'Ru101', 'Ru106',
       'Rh103', 'Pd105', 'Pd107', 'Ag109', 'Te132', 'I129', 'I131', 'Xe131',
       'Xe135', 'Cs133', 'Cs134', 'Cs135', 'Cs137', 'La139', 'Ce142', 'Nd143',
       'Nd145', 'Pm147', 'Sm149', 'Sm151']
all_nuc = actinides + fps

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

###############################################################################
#     Reaction Rates to direct tally for hybrid (Salcedo-Perez 2019 M&C)
###############################################################################
rr2 = ['fission', '(n,gamma)']
rr1 = ['(n,gamma)']

###############################################################################
#                  Initialize depletion calculation constants
###############################################################################
chain_file = '../data/depletion/chain_casl_pwr.xml'
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
op = openmc.deplete.CoupledOperator(model, chain_file=chain_file,
                                    fission_yield_mode="average",
                                    reaction_rate_mode='direct')
run(op)
mv_results('results/direct/depletion_results.h5')

###############################################################################
#                Run depletion with all nuclides flux tallied
###############################################################################
op = setup_flux_operator(reations=None, nuclides=None)
run(op)
mv_results('results/flux/depletion_results.h5')

###############################################################################
#                          Run hybrid depletion
###############################################################################
# pass
