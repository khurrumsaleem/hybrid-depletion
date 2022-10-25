from math import pi

import matplotlib.pyplot as plt
import openmc.deplete

###############################################################################
#                       List of Nuclides (Romano 2021)
###############################################################################
actinides = ['U234', 'U235', 'U236', 'U238', 'U239', 'Np239',
             'Pu238', 'Pu239', 'Pu240', 'Pu241', 'Pu242',
             'Am241', 'Am242', 'Am242_m1', 'Am243', 'Am244',
             'Cm242', 'Cm243', 'Cm244', 'Cm245', 'Cm246']
fps = ['Kr85', 'Sr90', 'Y90', 'Zr93', 'Mo95', 'Mo97', 'Tc99', 'Ru101', 'Ru106',
       'Rh103', 'Pd105', 'Pd107', 'Ag109', 'Te132', 'I129', 'I131', 'Xe131',
       'Xe135', 'Cs133', 'Cs134', 'Cs135', 'Cs137', 'La139', 'Ce142', 'Nd143',
       'Nd145', 'Pm147', 'Sm149', 'Sm151']
all_nuc = actinides + fps

radius = 0.39218
volume = pi * radius**2
day = 24*60*60

###############################################################################
#                              Load Results
###############################################################################
res_dir = openmc.deplete.Results("results/direct/depletion_results.h5")
res_flux = openmc.deplete.Results("results/flux/depletion_results.h5")

###############################################################################
#                      Plot Absolute Results per nuclide
###############################################################################
for nuc in all_nuc:
    # direct tally results
    _, atoms_dir = res_dir.get_atoms('1', nuc)
    conc_dir = atoms_dir * 1e-24/volume  # [atoms] [cm^2/b] / [cm^2] = atom/b

    # flux tally results
    time, atoms_flux = res_flux.get_atoms('1', nuc)
    conc_flux = atoms_flux * 1e-24/volume  # [atoms] [cm^2/b] / [cm^2] = atom/b

    fig, ax = plt.subplots()
    ax.plot(time/day, conc_dir, 'bo', label=f"{nuc} (direct)")
    ax.plot(time/day, conc_flux, 'kx', label=f"{nuc} (flux)")
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Atom/barn")
    ax.legend()
    ax.grid(True, which='both')

    plt.savefig("results/figures/abs_nuclides/{}.png".format(nuc))
    plt.close()
