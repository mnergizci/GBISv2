%% run_GBIS.m

addpath(genpath('D:/1.Phd_project/1.GBISv2/GBIS'));

testGBISPath
GBISrun('TR_sparse_exercise_empr.inp',[],'n','y','Z',1e6, 'n', 'n', 'n');
plot_grind('invert_ENU_Z','TR_sparse_exercise_empr/');
invResFile='TR_sparse_exercise_empr/invert_ENU_Z/';
generateFinalReport(invResFile);