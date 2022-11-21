load data/1_filtered.sdf
pseudoatom 1_mean, pos=[4.192599384597616, 64.08139510333896, 0.23697085474097038]
load data/2_filtered.sdf
pseudoatom 2_mean, pos=[20.331823215709246, 57.20518171750834, 3.66021271266139]
load data/3_filtered.sdf
pseudoatom 3_mean, pos=[16.756107854962135, 57.68595846349061, 2.4548977926088833]
load data/5_filtered.sdf
pseudoatom 5_mean, pos=[4.796297723562841, 61.83411852028851, -24.60506835344093]
load data/7_filtered.sdf
pseudoatom 7_mean, pos=[-15.29578158499439, 63.53265424015794, -12.02462511913339]
load data/8_filtered.sdf
pseudoatom 8_mean, pos=[29.80508338345017, 53.388053836281564, -3.598936685625911]

show spheres, 1_mean
show spheres, 2_mean
show spheres, 3_mean
show spheres, 5_mean
show spheres, 7_mean
show spheres, 8_mean

order *,yes
group glide, 1_mean 1_filtered 2_mean 2_filtered 3_mean 3_filtered 5_mean 5_filtered 7_mean 7_filtered 8_mean 8_filtered

load data/fda_eqina.sdf
load data/enamine_best_1000.sdf
load data/enamine_random_1000.sdf
load data/site_3_full_gnina_no_duplicates.sdf


set all_states, on
load data/3QE2_1_reduced.pdb
remove solvent