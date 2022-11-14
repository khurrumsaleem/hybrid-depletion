import openmc
from beavrs.builder import BEAVRS
import beavrs.constants as c

# model from here: https://github.com/mit-crpg/BEAVRS/blob/master/models/openmc/extract-assm.ipynb

b = BEAVRS()

# Get all OpenMC Lattices in a Python list
all_latts = b.openmc_geometry.get_all_lattices()
for id, latt in all_latts.items():
    if latt.name == "Fuel 1.6% enr instr no BAs":
        assembly = latt

# Create surface objects for our "root" cell"
lattice_sides = openmc.model.rectangular_prism(17*c.pinPitch, 17*c.pinPitch,
                                                   boundary_type='reflective')
min_z = openmc.ZPlane(z0=c.struct_LowestExtent, boundary_type='vacuum')
max_z = openmc.ZPlane(z0=c.struct_HighestExtent, boundary_type='vacuum')

# Create a "root" cell filled by the fuel assembly
root_cell = openmc.Cell(name='root cell',
                        fill=assembly,
                        region=lattice_sides & +min_z & -max_z
                       )

# Create a "root" universe with ID=0 and add the "root" cell to it
root_univ = openmc.Universe(name='root universe', cells=[root_cell])

# Create an OpenMC Geometry around root Universe
sub_geometry = openmc.Geometry(root_univ)

# Export the OpenMC Geometry to a "geometry.xml" file
sub_geometry.export_to_xml()

# Get a list of all OpenMC Materials
all_materials = sub_geometry.get_all_materials()

# Create a MaterialsFile object
materials = openmc.Materials(all_materials.values())

# Export the OpenMC Materials to a "materials.xml" file
materials.export_to_xml()

settings = openmc.Settings()

# Set any settings of interest
settings.batches = 80
settings.inactive = 20
settings.particles = int(5e5)

# set volume calc
mat_doms = list(sub_geometry.get_all_materials().values())
bb = sub_geometry.bounding_box
mat_calc = openmc.VolumeCalculation(domains=mat_doms,
                                    samples=100000000,
                                    lower_left=bb[0],
                                    upper_right=bb[1])
settings.volume_calculations = [mat_calc]

# Use a bounding box to define the starting source distribution
lower_left = [-17*c.pinPitch/2, -17*c.pinPitch/2, c.fuel_ActiveFuel_bot]
upper_right = [+17*c.pinPitch/2, +17*c.pinPitch/2, c.fuel_ActiveFuel_top]
settings.source = openmc.source.Source(
    openmc.stats.Box(lower_left, upper_right, only_fissionable=True))

settings.run_mode = 'volume'

# Export the settings to a "settings.xml" file
settings.export_to_xml()

# run volume calc
openmc.run()

# load volume info
vol_info = openmc.VolumeCalculation.from_hdf5('volume_1.h5')
for mat in materials:
    mat.add_volume_information(vol_info)

# re-export with volume information
materials.export_to_xml()