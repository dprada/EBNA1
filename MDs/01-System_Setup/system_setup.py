import warnings
warnings.filterwarnings("ignore")

import os
import molmodmt as m3t
import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as unit

# Loading 1B3T PDB file

system = m3t.load('pdb:1B3T','pdbfixer')

# Fixing PDB file

system = m3t.fix_pdb_structure(system)

# Removing crystal waters

system = m3t.remove_solvent(system)

# Adding temporary hydrogens and protonation states

system = m3t.add_hydrogens(system, pH=7.4, verbose=True)

# Building solvated truncated octahedral box
# The system is protonated again by LEaP

system = m3t.convert(system, 'openmm.Modeller')
system = m3t.solvate(system, add_hydrogens=True,
                     box_geometry="truncated_octahedral", clearance=14.0*unit.angstroms,
                     water='TIP3P', num_anions="neutralize", num_cations="neutralize",
                     to_form="openmm.Modeller", engine="LEaP")

# Fixing residue and atom names
# This section is unnecesary in next OpenMM versions
# See: https://github.com/pandegroup/openmm/pull/2314

for residue in system.topology.residues():
    if residue.name=='DG5':
        residue.name='DG'
    if residue.name=='DC3':
        residue.name='DC'

for atom in system.topology.atoms():
    if atom.name in ["H2'1", "H2'2", "H5'1", "H5'2"]:
        if atom.name=="H2'1":
            atom.name="H2'"
        elif atom.name=="H2'2":
            atom.name="H2''"
        elif atom.name=="H5'1":
            atom.name="H5'"
        elif atom.name=="H5'2":
            atom.name="H5''"

m3t.convert(system,'auxfile.pdb')
system = m3t.convert('auxfile.pdb', 'openmm.Modeller')
os.remove('auxfile.pdb')

# Output as PDB file

m3t.convert(system, 'system_init.pdb')

# Print out info
print('Num atoms')
print('Protonation state')
print('Volume box')
print('Num waters')
print('Num ions')

