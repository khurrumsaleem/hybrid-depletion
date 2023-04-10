from math import pi
import matplotlib.pyplot as plt
import openmc.deplete
import numpy as np


groups = [300, 500, 2500, 10000]
nuc_sets = ['all', 'mix', 'actinides']

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
_, keff_dir = direct_results.get_keff()

for egroup in groups:
    flux_results = openmc.deplete.Results(root_path.format(run_info=f"flux/{egroup}"))
    _, keff_flux = flux_results.get_keff()

    for nucs in nuc_sets:
        print(f"Plotting K-Effective {egroup} {nucs}")
        hybrid1_results = openmc.deplete.Results(root_path.format(run_info=f"hybrid1/{nucs}/{egroup}"))
        hybrid2_results = openmc.deplete.Results(root_path.format(run_info=f"hybrid2/{nucs}/{egroup}"))
        print("... loading complete")

        ###############################################################################
        #                      Plot K-eff regardless
        ###############################################################################

        _, keff_hy1 = hybrid1_results.get_keff()
        time, keff_hy2 = hybrid2_results.get_keff()

        fig, ax = plt.subplots()
        ax.errorbar(time/day, keff_dir[:,0], yerr=2 * abs(keff_dir[:,1]),
                    fmt='bo', ecolor='black', label='direct', capsize=3)
        ax.errorbar(time/day, keff_flux[:,0], yerr=2 * abs(keff_flux[:,1]),
                    fmt='kx', ecolor='black', label='flux', capsize=3)
        ax.errorbar(time/day, keff_hy1[:,0], yerr=2 * abs(keff_hy1[:,1]),
                    fmt='g+', ecolor='black', label='hybrid 1', capsize=3)
        ax.errorbar(time/day, keff_hy2[:,0], yerr=2 * abs(keff_hy2[:,1]),
                    fmt='m*', ecolor='black', label='hybrid 2', capsize=3)

        #ax.axhline(color='k', linestyle='--')
        ax.set_xlabel("Time [days]")
        ax.set_ylabel("k_eff")
        ax.set_title(f"K-effective, E={egroup}, {nucs}")
        ax.grid(True)
        ax.legend()
        plt.tight_layout()
        plt.savefig(f"figures/keff/keff_{nucs}_{egroup}.png")
        plt.close()

        h1_diff = (keff_hy1[:,0] - keff_dir[:,0]) * 1e5
        h2_diff = (keff_hy2[:,0] - keff_dir[:,0]) * 1e5
        h1_diff_err = np.sqrt((keff_hy1[:,1]**2 + keff_dir[:,1]**2)) * 1e5
        h2_diff_err = np.sqrt((keff_hy2[:,1]**2 + keff_dir[:,1]**2)) * 1e5
        fig, ax = plt.subplots()
        ax.errorbar(time/day, h1_diff, yerr=2 * abs(h1_diff_err),
                    fmt='g+', ecolor='black', label='hybrid 1', capsize=3)
        ax.errorbar(time/day, h2_diff, yerr=2 * abs(h2_diff_err),
                    fmt='m*', ecolor='black', label='hybrid 2', capsize=4)
        ax.set_xlabel("Time [days]")
        ax.set_ylabel("k_eff diff, pcm")
        ax.set_title(f"K-effective diff, E={egroup}, {nucs}")
        ax.legend()
        ax.grid(True, which='both')
        plt.tight_layout()
        plt.savefig(f"figures/keff/keff_{nucs}_{egroup}_diff.png")
        plt.close()

