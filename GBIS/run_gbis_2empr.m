%% run_GBIS.m

addpath(genpath('D:/1.Phd_project/1.GBISv2/GBIS'));

testGBISPath
GBISrun('TR_sparse_exercise_empr.inp',[],'n','y','Z',2e5, 'n', 'n', 'n');
plot_grind('invert_1_2_3_4_5_6_7_8_9_Z','/nobackup/mnergizci/new_GBIS/GBIS_V2.0.1/Example/sill_nugget_range/TR_sparse_exercise_empr/');
invResFile='/nobackup/mnergizci/new_GBIS/GBIS_V2.0.1/Example/sill_nugget_range/TR_sparse_exercise_empr/invert_1_2_3_4_5_6_7_8_9_Z/';
generateFinalReport(invResFile,100000);
generate_plots(invResFile,100000);
