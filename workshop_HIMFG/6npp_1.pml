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
     0.019482266,    0.996740818,   -0.078172363,\
     0.726869345,    0.039567064,    0.685619771,\
     0.686485469,   -0.070178606,   -0.723737776,\
     0.000202328,    0.000169608, -176.117889404,\
   -11.571586609,   -9.674448013,   -6.626377106,\
  -166.005874634,  518.257263184,   20.000000000 )

#show cartoon, protein_6NPP
#set cartoon_color, c_binding_domain, protein_6NPP
show surface, protein_6NPP
set surface_color, c_binding_domain, protein_6NPP
#set transparency, 0.8

show cartoon, dna_1B3T
color c_dna, dna_1B3T
set cartoon_transparency, 0.5, dna_1B3T

rebuild

## Quita la neblina
set depth_cue, 0

### Quito o pongo perspectiva
set orthoscopic, on

## Elijo últimas opciones edición para renderizar
set antialias, 2
set opaque_background, off
set ray_trace_gain, 0.1
set ray_trace_mode, 0
set ray_shadows, 1

# Renderizamos y guardamos
ray 1500,1200
png docking_1.png, dpi=300

## Para ejecutar el script desde la terminal:
## >>> pymol -qc monomero.pml

