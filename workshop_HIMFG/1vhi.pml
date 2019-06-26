# Loading pdbs

load 1VHI.pdb
remove solvent and 1VHI

load 1B3T.pdb
remove solvent and 1B3T

# Alignment
align 1VHI, 1B3T


# Selections
select protein_1VHI, resid 469-607 and 1VHI
select protein_1B3T, resid 461-607 and 1B3T
select dna_1B3T, resid 101-218 and 1B3T
deselect

# Colors definition

set_color c_binding_domain, [52, 101, 164]



# Representations

hide everything

bg_color white

set_view (\
     0.017421240,   -0.708301306,    0.705697656,\
    -0.045909986,    0.704489470,    0.708226800,\
    -0.998790801,   -0.044736903,   -0.020244779,\
     0.000069354,   -0.000017666, -176.120300293,\
    10.205041885,    4.652202606,   12.896821022,\
   138.858825684,  213.392517090,  -20.000000000 )

show cartoon, protein_1VHI
set cartoon_color, c_binding_domain, protein_1VHI
show surface, protein_1VHI
set surface_color, white, protein_1VHI
set transparency, 0.8
rebuild


## Quita la neblina
set depth_cue, 0

## Quito o pongo perspectiva
set orthoscopic, on

# Elijo últimas opciones de edición para renderizar
set antialias, 2
set opaque_background, off
set ray_trace_gain, 0.1
set ray_trace_mode, 0
set ray_shadows, 1

# Renderizamos y guardamos
ray 1500,1200
png 1vhi.png, dpi=300

# Para ejecutar el script desde la terminal:
# >>> pymol -qc monomero.pml

