%% run_GBIS.m

addpath(genpath('/nobackup/mnergizci/new_GBIS/GBIS_V2.0.1/'));

testGBISPath
invResFile='/nobackup/mnergizci/new_GBIS/GBIS_V2.0.1/Example/sill_nugget_range/TR_sparse_exercise_empr/ramp_removed/invert_1_2_3_4_5_6_7_8_9_Z/';
invResFile='/nobackup/mnergizci/new_GBIS/GBIS_V2.0.1/Example/sill_nugget_range/TR_sparse_exercise_empr/other_comp/ramp_removed/invert_1_2_3_4_5_6_7_8_9_Z/';
generate_plots(invResFile,100000);
