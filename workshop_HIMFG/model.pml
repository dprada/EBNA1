# Loading pdbs

load 1B3T.pdb
remove solvent and 1B3T

load 11262_2014_1101_MOESM7_ESM.pdb, EBNA1_model
remove solvent and EBNA1_model

# Selections
select model_mon_A, chain A and EBNA1_model
select model_mon_B, chain B and EBNA1_model
select model_protein, (chain A or chain B) and EBNA1_model
select model_dna, (chain C or chain D) and EBNA1_model

select dna_1B3T, resid 101-218 and 1B3T
deselect


# Alignment
#align 1VHI, 1B3T
#align 1B3T, EBNA1_model
align model_dna, dna_1B3T

# Colors definition

set_color c_binding_domain, [52, 101, 164]
set_color c_dna, [244, 172, 69]



# Representations

hide everything

bg_color white

set_view (\
     0.017421240,   -0.708301306,    0.705697656,\
    -0.045909986,    0.704489470,    0.708226800,\
    -0.998790801,   -0.044736903,   -0.020244779,\
     0.000066456,    0.000044197, -331.618988037,\
    28.261356354,  -13.285377502,   14.197230339,\
   294.369598389,  368.903381348,  -20.000000000 )

show cartoon, model_protein
set cartoon_color, c_binding_domain, model_protein
show surface, model_protein
set surface_color, white, model_protein
set transparency, 0.8

show cartoon, dna_1B3T
show surface, dna_1B3T
color c_dna, dna_1B3T
set surface_color, white, dna_1B3T
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
png model.png, dpi=300

# Para ejecutar el script desde la terminal:
# >>> pymol -qc monomero.pml

