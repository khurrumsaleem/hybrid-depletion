import openmc
import openmc.deplete
import numpy as np

openmc.deplete.pool.NUM_PROCESSES = 1
openmc.deplete.pool.USE_MULTIPROCESSING = False

model = openmc.model.Model.from_xml()

###############################################################################
#                   Initialize and run depletion calculation
###############################################################################

# list of nuclides taken from Romano 2021
actinides = ['U234', 'U235', 'U236', 'U238', 'U239','Np239',
             'Pu238', 'Pu239', 'Pu240', 'Pu241', 'Pu242',
             'Am241', 'Am242', 'Am242_m1', 'Am243', 'Am244',
             'Cm242', 'Cm243', 'Cm244', 'Cm245', 'Cm246']
fps = []  # fission products

# groups structure from Salcedo-Perez 2019 M&C
groups500 = list(np.logspace(np.log10(1e-5), np.log10(400e3), 50, endpoint=False)) +\
    list(np.logspace(np.log10(400e3), np.log10(3.21e6), 100, endpoint=False)) +\
            list(np.logspace(np.log10(3.21e6), np.log10(8.025e6), 260)) +\
            list(np.logspace(np.log10(8.05e6), np.log10(20e6), 90))

# setup helper - I don't understand the reaction rates and nuclides thing
helper = openmc.deplete.helpers.FluxCollapseHelper(162, 2, groups500)

# Set up depletion operator
chain_file = '../data/depletion/chain_casl_pwr.xml'
op = openmc.deplete.CoupledOperator(model, chain_file=chain_file,
                                    fission_yield_mode="average",
                                    reaction_rate_mode='flux',
                                    reaction_rate_opts={'energies': groups500})

# cumulative steps in MWd/kg
burnup_cum = np.array([
    0.1, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0 #,
    #12.5, 15.0, 17.5, 20.0, 22.5, 25.0, 27.5, 30.0, 32.5, 35.0, 37.5,
    #40.0, 42.5, 45.0, 47.5, 50.0
])
burnup = np.diff(burnup_cum, prepend=0.0)
power = 174  # W/cm

# Perform simulation using the predictor algorithm
integrator = openmc.deplete.PredictorIntegrator(op, burnup, power,
                                                timestep_units='MWd/kg')
integrator.integrate()
