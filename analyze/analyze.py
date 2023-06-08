from math import pi
import matplotlib.pyplot as plt
import openmc.deplete
import argparse

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
all_nuclides = list(set(actinides + mix))  # plot all nuclides of interest regardless
if args.nuclides == 'all':
    nuc_set = list(set(actinides + mix))
elif args.nuclides == 'actinides':
    nuc_set = actinides
elif args.nuclides == 'mix':
    nuc_set = mix
egroup = args.groups

radius = 0.39218
volume = pi * radius**2
day = 24*60*60
vol_mult = 1e-24/volume

###############################################################################
#                              Load Results
###############################################################################
print("Loading results...")
root_path = "/home/kkiesling/depletion/hybrid-depletion/deplete/results/pin/comp-redo/{run_info}/depletion_results.h5"
direct_results = openmc.deplete.Results(root_path.format(run_info="direct"))
#flux_results = openmc.deplete.Results(root_path.format(run_info=f"flux/{egroup}"))
#hybrid1_results = openmc.deplete.Results(root_path.format(run_info=f"hybrid1/{args.nuclides}/{egroup}"))
hybrid2_results = openmc.deplete.Results(root_path.format(run_info=f"hybrid2/{args.nuclides}/{egroup}"))
print("... loading complete")

###############################################################################
#                Plot Absolute Results per nuclide compared to direct
###############################################################################
print("calculating nuclide diffs and plotting")
for nuc in sorted(all_nuclides):
    ast = "*" if nuc in nuc_set else "" # label if it was direct tallied for this hybrid data
    print(f"\tPlotting {nuc}{ast} results")

    # direct results
    _, atoms_dir = direct_results.get_atoms('1', nuc)
    conc_dir = atoms_dir * vol_mult # [atoms] [cm^2/b] / [cm^2] = atom/b

    # flux results
    #_, atoms_flux = flux_results.get_atoms('1', nuc)
    #conc_flux = atoms_flux * vol_mult  # [atoms] [cm^2/b] / [cm^2] = atom/b

    # hybrid 1  results
    #_, atoms_hy1 = hybrid1_results.get_atoms('1', nuc)
    #conc_hy1 = atoms_hy1 * vol_mult  # [atoms] [cm^2/b] / [cm^2] = atom/b

    # hybrid 2 results
    time, atoms_hy2 = hybrid2_results.get_atoms('1', nuc)
    conc_hy2 = atoms_hy2 * vol_mult  # [atoms] [cm^2/b] / [cm^2] = atom/b

    # plot absolute
    fig, ax = plt.subplots()
    ax.plot(time/day, conc_dir, 'bo', label="direct")
    #ax.plot(time/day, conc_flux, 'kx', label="flux")
    #ax.plot(time/day, conc_hy1, 'g+', label="hybrid 1")
    ax.plot(time/day, conc_hy2, 'm*', label="hybrid 2")
    ax.set_title(f"{nuc}{ast}, E={egroup}, {args.nuclides} - original")
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Concentration [atoms/b]")
    ax.legend()
    ax.grid(True, which='both')
    plt.tight_layout()
    plt.savefig(f"figures/{args.nuclides}/{egroup}/{nuc}-absolute-hy2.png")
    plt.close()

    # plot diff compared to direct
    #flux_diff = (conc_flux - conc_dir) #/ conc_dir * 100
    h2_diff = (conc_hy2 - conc_dir)    #/ conc_dir * 100
    #h1_diff = (conc_hy1 - conc_dir)    #/ conc_dir * 100

    # plot diff
    fig, ax = plt.subplots()
    #ax.plot(time/day, flux_diff, 'kx', label="flux")
    #ax.plot(time/day, h1_diff, 'g+', label="hybrid 1")
    ax.plot(time/day, h2_diff, 'm*', label="hybrid 2")
    ax.set_title(f"{nuc}{ast}, E={egroup}, {args.nuclides} - original")
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Concentration Diff [conc - direct]")
    ax.legend()
    ax.grid(True, which='both')
    plt.tight_layout()
    plt.savefig(f"figures/{args.nuclides}/{egroup}/{nuc}-diff-hy2.png")
    plt.close()

    # plot diff compared to direct
    #flux_diff = (conc_flux - conc_dir) / conc_dir * 100
    h2_diff = (conc_hy2 - conc_dir) / conc_dir * 100
    #h1_diff = (conc_hy1 - conc_dir) / conc_dir * 100

    # plot diff
    fig, ax = plt.subplots()
    #ax.plot(time/day, flux_diff, 'kx', label="flux")
    #ax.plot(time/day, h1_diff, 'g+', label="hybrid 1")
    ax.plot(time/day, h2_diff, 'm*', label="hybrid 2")
    ax.set_title(f"{nuc}{ast}, E={egroup}, {args.nuclides} - original")
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Relative Concentration [(conc - direct)/direct]")
    ax.legend()
    ax.grid(True, which='both')
    plt.tight_layout()
    plt.savefig(f"figures/{args.nuclides}/{egroup}/{nuc}-relative-hy2.png")
    plt.close()


# ###############################################################################
# #             Plot EOL concentration diffs compared to direct tally
# ###############################################################################
# def plot_nuc_diffs(nuclides, results, nuc_name, res_name):
#     eol = []
#     for nuc in sorted(nuclides):
#         # direct tally results
#         _, atoms_dir = res_dir.get_atoms('1', nuc)
#         conc_dir = atoms_dir * 1e-24/volume  # [atoms] [cm^2/b] / [cm^2] = atom/b

#         # flux tally results
#         _, atoms = results.get_atoms('1', nuc)
#         conc = atoms * 1e-24/volume  # [atoms] [cm^2/b] / [cm^2] = atom/b

#         conc_diff = (conc[-1] - conc_dir[-1]) #/ conc_dir[-1] * 100
#         eol.append(conc_diff)

#     fig, ax = plt.subplots(figsize=(5,10))
#     y_pos = np.arange(len(nuclides))
#     ax.barh(y_pos, eol, align='center')
#     ax.set_yticks(y_pos, labels=sorted(nuclides))
#     ax.axvline(color='k', linestyle='--')
#     ax.grid(True, axis='x')
#     ax.invert_yaxis()  # labels read top-to-bottom
#     ax.set_xlabel('EOL Difference')
#     ax.set_title(f"{res_name}, {nuc_name}")
#     ax.grid(True, which='both')
#     plt.tight_layout()
#     plt.savefig("figures/all/{}/eol_absolute_{}_{}.png".format(egroup, nuc_name.strip(" ").lower(), res_name.strip(" ").lower()))
#     plt.close()

# # plot_nuc_diffs(actinides, res_flux, 'Actinides', 'Flux')
# # plot_nuc_diffs(fps, res_flux, 'Fission Products', 'Flux')
# plot_nuc_diffs(nuclides, res_hy1, args.nuclides.capitalize(), 'Hybrid 1')
# # plot_nuc_diffs(fps, res_hy1, 'Fission Products', 'Hybrid 1')
# # plot_nuc_diffs(actinides, res_hy2, 'Actinides', 'Hybrid 2')
# # plot_nuc_diffs(fps, res_hy2, 'Fission Products', 'Hybrid 2')
