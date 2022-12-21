from math import pi
import matplotlib.pyplot as plt
import openmc.deplete
from uncertainties import unumpy as unp
import numpy as np
import argparse
import pathlib

###############################################################################
#                           Parse Args
###############################################################################
parser = argparse.ArgumentParser(description='Specify Simulations to run')
parser.add_argument('-g', '--groups', default=500, type=int, choices=[300, 500, 2500, 10000],
                    help='Number of energy groups to use; 500 or 300 or 2500 (default=500)')
parser.add_argument('-n', '--nuclides', default='all', type=str, choices=['all', 'actinides', 'mix'],
                    help="Which nuclides to direct tally")
args = parser.parse_args()

###############################################################################
#                       List of Nuclides
###############################################################################
actinides = ['Th230', 'Th231', 'Th232', 'Th234', 'Pa231', 'Pa232',
             'Pa233', 'U232', 'U233', 'U234', 'U235', 'U236', 'U237',
             'U238', 'U239', 'Np236', 'Np237', 'Np238', 'Np239',
             'Pu236', 'Pu238', 'Pu239', 'Pu240', 'Pu241', 'Pu242',
             'Pu243', 'Am241', 'Am242', 'Am242_m1', 'Am243', 'Cm242', 'Cm244']
mix = ['Kr83', 'Tc99', 'Ru101', 'Rh103', 'Rh105', 'Xe131', 'Xe133',
       'Xe135', 'Cs133', 'Cs135', 'Pr143', 'Nd143', 'Nd145',
       'Pm147', 'Pm148_m1', 'Pm149', 'Pm151', 'Sm149', 'Sm150',
       'Sm151', 'Sm152', 'Sm153', 'Eu153', 'Eu155', 'Gd157', 'U234',
       'U235', 'U236', 'U237', 'U238', 'Pu239', 'Pu240']
if args.nuclides == 'all':
    nuclides = list(set(actinides + mix))
elif args.nuclides == 'actinides':
    nuclides = actinides
elif args.nuclides == 'mix':
    nuclides = mix
egroup = args.groups

radius = 0.39218
volume = pi * radius**2
day = 24*60*60

###############################################################################
#                              Load Results
###############################################################################
print("loading results")
res_dir = openmc.deplete.Results("../deplete/results/pin/comprehensive/direct/300/depletion_results.h5")
res_flux = openmc.deplete.Results("../deplete/results/pin/comprehensive/flux/{}/depletion_results.h5".format(egroup))
res_hy1 = openmc.deplete.Results("../deplete/results/pin/comprehensive/hybrid1/{}/{}/depletion_results.h5".format(args.nuclides, egroup))
res_hy2 = openmc.deplete.Results("../deplete/results/pin/comprehensive/hybrid2/{}/{}/depletion_results.h5".format(args.nuclides, egroup))

###############################################################################
#                Plot Absolute Results per nuclide compared to direct
###############################################################################
print("calculating nuclide diffs and plotting")
for nuc in sorted(nuclides):
    # direct tally results
    _, atoms_dir = res_dir.get_atoms('1', nuc)
    conc_dir = atoms_dir * 1e-24/volume  # [atoms] [cm^2/b] / [cm^2] = atom/b

    # flux tally results
    #_, atoms_flux = res_flux.get_atoms('1', nuc)
    #conc_flux = atoms_flux * 1e-24/volume  # [atoms] [cm^2/b] / [cm^2] = atom/b

    # flux tally results
    time, atoms_hy1 = res_hy1.get_atoms('1', nuc)
    conc_hy1 = atoms_hy1 * 1e-24/volume  # [atoms] [cm^2/b] / [cm^2] = atom/b

    # flux tally results
    time2, atoms_hy2 = res_hy2.get_atoms('1', nuc)
    conc_hy2 = atoms_hy2 * 1e-24/volume  # [atoms] [cm^2/b] / [cm^2] = atom/b

    # # plot absolute
    # fig, ax = plt.subplots()
    # ax.plot(time/day, conc_dir, 'bo', label="direct")
    # ax.plot(time/day, conc_flux, 'kx', label="flux")
    # ax.plot(time/day, conc_hy1, 'g+', label="hybrid 1")
    # ax.plot(time/day, conc_hy2, 'm*', label="hybrid 2")
    # ax.set_title(f"{nuc}")
    # ax.set_xlabel("Time (days)")
    # ax.set_ylabel("Atom/barn")
    # ax.legend()
    # ax.grid(True, which='both')
    # plt.tight_layout()
    # plt.savefig("figures/{}/{}/{}.png".format(args.nuclides, egroup, nuc))
    # plt.close()

    # plot diff compared to direct
    #h1_diff = (conc_hy1 - conc_dir) #/ conc_dir * 100
    #flux_diff = (conc_flux - conc_dir) #/ conc_dir * 100
    #h2_diff = (conc_hy2 - conc_dir[:len(conc_hy2)]) #/ conc_dir * 100
    #print(nuc, h1_diff[-1], atoms_hy1[-1], atoms_dir[-1])
    #h2_diff = (atoms_hy2 - atoms_dir) #/ atoms_dir * 100

    fig, ax = plt.subplots()
    ax.plot(time/day, conc_dir, 'rx', label="direct")
    ax.plot(time/day, conc_hy1, 'g+', label="hybrid 1")
    ax.plot(time2/day, conc_hy2, 'm*', label="hybrid 2")
    #ax.plot(time/day, conc_flux, 'b.', label="flux")

    #ax.plot(time/day, h2_diff, 'm*', label="hybrid 2")
   # ax.axhline(color='k', linestyle='--')
    ax.set_title(f"{nuc}")
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Concentration [atoms/b]")
    ax.legend()
    ax.grid(True, which='both')
    plt.tight_layout()
    plt.savefig("figures/{}/{}/{}.png".format(args.nuclides, egroup, nuc))
    plt.close()

###############################################################################
#                      Plot K-eff
###############################################################################
def plot_keff_diff(results, name):
    # direct tally results
    _, keff_dir = res_dir.get_keff()
    # results to compare
    time, keff_hy1 = res_hy1.get_keff()
    time2, keff_hy2 = res_hy2.get_keff()
    time, keff_flux = res_flux.get_keff()

    #k_diff = (keff_hy1 - keff_dir) * 1e5
    fig, ax = plt.subplots()
    ax.errorbar(time/day, keff_dir[:,0], yerr=2 * abs(keff_dir[:,1]),
                fmt='rx', ecolor='black', label='direct', markersize=8, capsize=3)
    ax.errorbar(time/day, keff_hy1[:,0], yerr=2 * abs(keff_hy1[:,1]),
                fmt='g+', ecolor='black', label='hybrid 1', markersize=8, capsize=3)
    ax.errorbar(time2/day, keff_hy2[:,0], yerr=2 * abs(keff_hy2[:,1]),
                fmt='m*', ecolor='black', label='hybrid 2', markersize=8, capsize=3)
    ax.errorbar(time/day, keff_flux[:,0], yerr=2 * abs(keff_flux[:,1]),
                fmt='b.', ecolor='black', label='flux', markersize=8, capsize=3)
    #ax.axhline(color='k', linestyle='--')
    ax.set_xlabel("Time [days]")
    ax.set_ylabel("k_eff")
    ax.grid(True)
    ax.legend()
    plt.tight_layout()
    plt.savefig("figures/keff_{}.png".format(egroup))
    plt.close()

# plot_keff_diff(res_flux, 'flux')
plot_keff_diff(res_hy1, 'hybrid 1')
# plot_keff_diff(res_hy2, 'hybrid 2')

###############################################################################
#             Plot EOL concentration diffs compared to direct tally
###############################################################################
def plot_nuc_diffs(nuclides, results, nuc_name, res_name):
    eol = []
    for nuc in sorted(nuclides):
        # direct tally results
        _, atoms_dir = res_dir.get_atoms('1', nuc)
        conc_dir = atoms_dir * 1e-24/volume  # [atoms] [cm^2/b] / [cm^2] = atom/b

        # flux tally results
        _, atoms = results.get_atoms('1', nuc)
        conc = atoms * 1e-24/volume  # [atoms] [cm^2/b] / [cm^2] = atom/b

        conc_diff = (conc[-1] - conc_dir[-1]) #/ conc_dir[-1] * 100
        eol.append(conc_diff)

    fig, ax = plt.subplots(figsize=(5,10))
    y_pos = np.arange(len(nuclides))
    ax.barh(y_pos, eol, align='center')
    ax.set_yticks(y_pos, labels=sorted(nuclides))
    ax.axvline(color='k', linestyle='--')
    ax.grid(True, axis='x')
    ax.invert_yaxis()  # labels read top-to-bottom
    ax.set_xlabel('EOL Difference')
    ax.set_title(f"{res_name}, {nuc_name}")
    ax.grid(True, which='both')
    plt.tight_layout()
    plt.savefig("figures/all/{}/eol_absolute_{}_{}.png".format(egroup, nuc_name.strip(" ").lower(), res_name.strip(" ").lower()))
    plt.close()

# plot_nuc_diffs(actinides, res_flux, 'Actinides', 'Flux')
# plot_nuc_diffs(fps, res_flux, 'Fission Products', 'Flux')
plot_nuc_diffs(nuclides, res_hy1, args.nuclides.capitalize(), 'Hybrid 1')
# plot_nuc_diffs(fps, res_hy1, 'Fission Products', 'Hybrid 1')
# plot_nuc_diffs(actinides, res_hy2, 'Actinides', 'Hybrid 2')
# plot_nuc_diffs(fps, res_hy2, 'Fission Products', 'Hybrid 2')
