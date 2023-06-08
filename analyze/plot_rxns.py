from math import pi
import matplotlib.pyplot as plt
import openmc.deplete
import argparse


radius = 0.39218
volume = pi * radius**2
day = 24*60*60
vol_mult = 1e-24/volume

parser = argparse.ArgumentParser(description='Specify nuclides and rxns')
parser.add_argument('nuclide', type=str, help="Which nuclide to plot")
parser.add_argument('rxn', type=str, help="reaction rate to plot")

args = parser.parse_args()

nuc = args.nuclide
rxn = args.rxn
rxn_name = rxn.strip('(').strip(')').replace(',', '_')
###############################################################################
#                              Load Results
###############################################################################
print("Loading results...")
root_path = "/home/kkiesling/depletion/hybrid-depletion/deplete/results/pin/comp-redo/{run_info}/depletion_results.h5"
direct_results = openmc.deplete.Results(root_path.format(run_info="direct"))
hybrid2_mod_results = openmc.deplete.Results(root_path.format(run_info=f"hybrid2-mod/actinides/500"))
hybrid2_results = openmc.deplete.Results(root_path.format(run_info=f"hybrid2/actinides/500"))
print("... loading complete")

# get reaction rates
time, dire = direct_results.get_reaction_rate("1", nuc, rxn)
_, hy2 = hybrid2_results.get_reaction_rate("1", nuc, rxn)
_, hy2_mod = hybrid2_mod_results.get_reaction_rate("1", nuc, rxn)


# plot absolute
fig, ax = plt.subplots()
ax.plot(time/day, dire, 'bo', label="direct")
ax.plot(time/day, hy2, 'm*', label="hybrid 2 - original")
ax.plot(time/day, hy2_mod, 'g+', label="hybrid 2 - modified")
ax.set_title(f"{nuc}, {rxn}")
ax.set_xlabel("Time (days)")
ax.set_ylabel("Reaction Rate, absolute")
ax.legend()
ax.grid(True, which='both')
plt.tight_layout()
plt.savefig(f"figures/actinides/500/rxns/{nuc}-{rxn_name}-abs.png")
plt.close()

hy2_rel = (hy2 - dire) / dire * 100
hy2_mod_rel = (hy2_mod - dire) / dire * 100

fig, ax = plt.subplots()
#ax.plot(time/day, hy2_rel, 'm*', label="hybrid 2 - original")
ax.plot(time/day, hy2_mod_rel, 'g+', label="hybrid 2 - modified")
ax.set_title(f"{nuc}, {rxn}, relative")
ax.set_xlabel("Time (days)")
ax.set_ylabel("Reaction Rate, relative difference [%]")
ax.legend()
ax.grid(True, which='both')
plt.tight_layout()
plt.savefig(f"figures/actinides/500/rxns/{nuc}-{rxn_name}-rel-modonly.png")
plt.close()

