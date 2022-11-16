from math import pi
import matplotlib.pyplot as plt
import openmc.deplete
from uncertainties import unumpy as unp
import numpy as np

###############################################################################
#                       List of Nuclides + constants (Romano 2021)
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
res_hy1 = openmc.deplete.Results("../deplete/results/hybrid1/depletion_results.h5")
res_hy2 = openmc.deplete.Results("../deplete/results/hybrid2/depletion_results.h5")

###############################################################################
#                Plot Absolute Results per nuclide compared to direct
###############################################################################
for nuc in all_nuc:
    # direct tally results
    _, atoms_dir = res_dir.get_atoms('1', nuc)
    conc_dir = atoms_dir * 1e-24/volume  # [atoms] [cm^2/b] / [cm^2] = atom/b

    # flux tally results
    _, atoms_flux = res_flux.get_atoms('1', nuc)
    conc_flux = atoms_flux * 1e-24/volume  # [atoms] [cm^2/b] / [cm^2] = atom/b

    # flux tally results
    _, atoms_hy1 = res_hy1.get_atoms('1', nuc)
    conc_hy1 = atoms_hy1 * 1e-24/volume  # [atoms] [cm^2/b] / [cm^2] = atom/b

    # flux tally results
    time, atoms_hy2 = res_hy2.get_atoms('1', nuc)
    conc_hy2 = atoms_hy2 * 1e-24/volume  # [atoms] [cm^2/b] / [cm^2] = atom/b

    # plot absolute
    fig, ax = plt.subplots()
    ax.plot(time/day, conc_dir, 'bo', label="direct")
    ax.plot(time/day, conc_flux, 'kx', label="flux")
    ax.plot(time/day, conc_hy1, 'g+', label="hybrid 1")
    ax.plot(time/day, conc_hy2, 'm*', label="hybrid 2")
    ax.set_title(f"{nuc}")
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Atom/barn")
    ax.legend()
    ax.grid(True, which='both')
    plt.tight_layout()
    plt.savefig("figures/abs_nuclides/{}.png".format(nuc))
    plt.close()

    # plot diff compared to direct
    h1_diff = (atoms_hy1 - atoms_dir) / atoms_dir * 100
    h2_diff = (atoms_hy2 - atoms_dir) / atoms_dir * 100

    fig, ax = plt.subplots()
    ax.plot(time/day, h1_diff, 'g+', label="hybrid 1")
    ax.plot(time/day, h2_diff, 'm*', label="hybrid 2")
    ax.axhline(color='k', linestyle='--')
    ax.set_title(f"{nuc}")
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Difference compared to direct (%)")
    ax.legend()
    ax.grid(True, which='both')
    plt.tight_layout()
    plt.savefig("figures/diff_nuclides/{}.png".format(nuc))
    plt.close()

###############################################################################
#                      Plot K-eff
###############################################################################
def plot_keff_diff(results, name):
    # direct tally results
    _, keff_dir = res_dir.get_keff()
    # results to compare
    time, keff = results.get_keff()

    k_diff = (keff - keff_dir) * 1e5
    fig, ax = plt.subplots()
    ax.errorbar(time/day, k_diff[:,0], yerr=2 * abs(k_diff[:,1]),
                fmt='b.', ecolor='black')
    ax.axhline(color='k', linestyle='--')
    ax.set_xlabel("Time [days]")
    ax.set_ylabel("k_{} - k_direct [pcm]".format(name))
    ax.grid(True)
    plt.savefig("figures/keff_diff_{}.png".format(name.strip(" ").lower()))
    plt.close()

plot_keff_diff(res_flux, 'flux')
plot_keff_diff(res_hy1, 'hybrid 1')
plot_keff_diff(res_hy2, 'hybrid 2')

###############################################################################
#             Plot EOL concentration diffs compared to direct tally
###############################################################################
def plot_nuc_diffs(nuclides, results, nuc_name, res_name):
    eol = []
    for nuc in nuclides:
        # direct tally results
        _, atoms_dir = res_dir.get_atoms('1', nuc)
        conc_dir = atoms_dir * 1e-24/volume  # [atoms] [cm^2/b] / [cm^2] = atom/b

        # flux tally results
        _, atoms = results.get_atoms('1', nuc)
        conc = atoms * 1e-24/volume  # [atoms] [cm^2/b] / [cm^2] = atom/b

        conc_diff = (conc - conc_dir) / conc_dir * 100
        eol.append(conc_diff[-1])

    fig, ax = plt.subplots(figsize=(5,8))
    y_pos = np.arange(len(nuclides))
    ax.barh(y_pos, eol, align='center')
    ax.set_yticks(y_pos, labels=nuclides)
    ax.axvline(color='k', linestyle='--')
    ax.grid(True, axis='x')
    ax.invert_yaxis()  # labels read top-to-bottom
    ax.set_xlabel('EOL Difference (%)')
    ax.set_title(f"{res_name}, {nuc_name}")
    plt.tight_layout()
    plt.savefig("figures/eol_{}_{}.png".format(nuc_name.strip(" ").lower(), res_name.strip(" ").lower()))
    plt.close()

plot_nuc_diffs(actinides, res_flux, 'Actinides', 'Flux')
plot_nuc_diffs(fps, res_flux, 'Fission Products', 'Flux')
plot_nuc_diffs(actinides, res_hy1, 'Actinides', 'Hybrid 1')
plot_nuc_diffs(fps, res_hy1, 'Fission Products', 'Hybrid 1')
plot_nuc_diffs(actinides, res_hy2, 'Actinides', 'Hybrid 2')
plot_nuc_diffs(fps, res_hy2, 'Fission Products', 'Hybrid 2')
