# Loading pdbs

#load 1VHI.pdb
#remove solvent and 1VHI

load 1B3T.pdb
remove solvent and 1B3T

load 6npp.pdb, 6NPP
remove solvent and 6NPP


# Alignment
align 1B3T, 6NPP


# Selections
select protein_6NPP, resid 467-607 and 6NPP
select ligand_6NPP, resn KWG

#select dna_1B3T, resname KWG
deselect

## Colors definition

set_color c_binding_domain, [52, 101, 164]
set_color c_dna, [244, 172, 69]



## Representations

hide everything

bg_color white

set_view (\
     0.017421240,   -0.708301306,    0.705697656,\
    -0.045909986,    0.704489470,    0.708226800,\
    -0.998790801,   -0.044736903,   -0.020244779,\
     0.000069354,   -0.000017666, -176.120300293,\
    10.205041885,    4.652202606,   12.896821022,\
   138.858825684,  213.392517090,  -20.000000000 )

show cartoon, protein_6NPP
set cartoon_color, c_binding_domain, protein_6NPP
show surface, protein_6NPP
set surface_color, white, protein_6NPP
set transparency, 0.8

#show cartoon, dna_1B3T
#show surface, dna_1B3T
#color c_dna, dna_1B3T
#set surface_color, white, dna_1B3T
#set transparency, 0.8

rebuild

## Quita la neblina
set depth_cue, 0

### Quito o pongo perspectiva
set orthoscopic, on

## Elijo últimas opciones de edición para renderizar
#set antialias, 2
#set opaque_background, off
#set ray_trace_gain, 0.1
#set ray_trace_mode, 0
#set ray_shadows, 1
#
## Renderizamos y guardamos
#ray 1500,1200
#png 1b3t.png, dpi=300
#
## Para ejecutar el script desde la terminal:
## >>> pymol -qc monomero.pml
#
