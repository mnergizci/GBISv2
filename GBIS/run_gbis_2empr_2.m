%% run_GBIS.m

addpath(genpath('/nobackup/mnergizci/new_GBIS/GBIS_V2.0.1/'));

testGBISPath

GBISrun('TR_sparse_exercise_empr.inp',[1,2,4,5,7,8],'n','Z',8e5, 'n', 'n', 'y');
plot_grind('invert_1_2_4_5_7_8_Z','/nobackup/mnergizci/new_GBIS/GBIS_V2.0.1/Example/sill_nugget_range/TR_sparse_exercise_empr/');
invResFile='/nobackup/mnergizci/new_GBIS/GBIS_V2.0.1/Example/sill_nugget_range/TR_sparse_exercise_empr/invert_1_2_4_5_7_8_Z'
generateFinalReport(invResFile,100000);
generate_plots(invResFile,100000);
