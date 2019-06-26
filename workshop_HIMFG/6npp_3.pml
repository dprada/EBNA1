# Loading pdbs

load 1B3T.pdb
remove solvent and 1B3T

load 6npp.pdb, 6NPP
remove solvent and 6NPP

# Alignment
align 1B3T, 6NPP

# Selections
select protein_6NPP, resid 467-607 and 6NPP
select ligand_6NPP, resn KWG
select protein_1B3T, resid 461-607 and 1B3T
select dna_1B3T, resid 101-218 and 1B3T

deselect

## Colors definition

set_color c_binding_domain, [52, 101, 164]
set_color c_dna, [244, 172, 69]

## Representations

hide everything

bg_color white

set_view (\
     0.125951633,    0.982413232,   -0.137722805,\
     0.946368158,   -0.160621241,   -0.280289769,\
    -0.297487885,   -0.095035166,   -0.949971437,\
     0.000504859,    0.000122059, -165.274459839,\
   -11.373695374,  -24.453845978,   -2.146590710,\
   136.197769165,  194.336532593,   20.000000000 )

show cartoon, protein_6NPP
set cartoon_color, c_binding_domain, protein_6NPP
show surface, protein_6NPP
set surface_color, red, protein_6NPP
set transparency, 0.4
set surface_cavity_mode, 1

show cartoon, dna_1B3T
color c_dna, dna_1B3T

rebuild

## Quita la neblina
#set depth_cue, 0

### Quito o pongo perspectiva
set orthoscopic, on

### Elijo últimas opciones edición para renderizar
set antialias, 2
set opaque_background, off
set ray_trace_gain, 0.1
set ray_trace_mode, 0
set ray_shadows, 1

## Renderizamos y guardamos
ray 1500,1200
png docking_3.png, dpi=300

## Para ejecutar el script desde la terminal:
## >>> pymol -qc monomero.pml

