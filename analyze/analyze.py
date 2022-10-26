from math import pi
import matplotlib.pyplot as plt
import openmc.deplete
from uncertainties import unumpy as unp
import numpy as np

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
res_dir = openmc.deplete.Results("../deplete/results/direct/depletion_results.h5")
res_flux = openmc.deplete.Results("../deplete/results/flux/depletion_results.h5")

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

    plt.savefig("figures/abs_nuclides/{}.png".format(nuc))
    plt.close()

    fig, ax = plt.subplots()
    ax.plot(time/day, (conc_flux - conc_dir)/conc_dir, 'bo', label=f"{nuc}")
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Relative Difference")
    ax.legend()
    ax.grid(True, which='both')

    plt.savefig("figures/comp_nuclides/{}.png".format(nuc))
    plt.close()

###############################################################################
#                      Plot K-eff
###############################################################################
# direct tally results
_, keff_dir = res_dir.get_keff()
# flux tally results
time, keff_flux = res_flux.get_keff()

k_flux = unp.uarray(keff_flux[:, 0], keff_flux[:, 1])
k_dir = unp.uarray(keff_dir[:, 0], keff_dir[:, 1])
k_diff = (keff_flux - keff_dir) * 1e5
fig, ax = plt.subplots()
ax.errorbar(time/day, k_diff[:,0], yerr=2 * abs(k_diff[:,1]),
            fmt='b.', ecolor='black')
ax.axhline(color='k', linestyle='--')
ax.set_xlabel("Time [days]")
ax.set_ylabel("$k_{flux} - k_{direct}$ [pcm]")
ax.grid(True)
plt.savefig("figures/keff_diff.png")
plt.close()

###############################################################################
#                      Plot EOL and BOL concentration diffs
###############################################################################

def plot_nuc_diffs(nuclides, name):
    eol = []
    for nuc in nuclides:
        # direct tally results
        _, atoms_dir = res_dir.get_atoms('1', nuc)
        conc_dir = atoms_dir * 1e-24/volume  # [atoms] [cm^2/b] / [cm^2] = atom/b

        # flux tally results
        time, atoms_flux = res_flux.get_atoms('1', nuc)
        conc_flux = atoms_flux * 1e-24/volume  # [atoms] [cm^2/b] / [cm^2] = atom/b

        conc_diff = (conc_flux - conc_dir) / conc_dir
        eol.append(conc_diff[-1])

    fig, ax = plt.subplots()
    y_pos = np.arange(len(nuclides))
    ax.barh(y_pos, eol, align='center')
    ax.set_yticks(y_pos, labels=nuclides)
    #ax.invert_yaxis()  # labels read top-to-bottom
    ax.set_xlabel('EOL Difference')
    ax.set_title(f"{name}")
    plt.savefig("figures/eol_{}.png".format(name.strip(" ").lower()))
    plt.show()

plot_nuc_diffs(actinides, 'Actinides')
plot_nuc_diffs(fps, 'Fission Products')
