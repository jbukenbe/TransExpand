function array_call(array_id)

%% Load Data
load tstnep_param_data param_set
tstnep_params = param_set{array_id};


%% Run Two Stage TNEP problem
TSTNEP_problem = OptTSTNEPClass(tstnep_params);
[TSTNEP_problem, opt_out] = TSTNEP_problem.solve();


%% Save Data
outfile_name = sprintf('%s_%d','tstnep_run',array_id);
m = matfile(outfile_name);

m.tstnep_problem = TSTNEP_problem; 
m.opt_out = opt_out;
m.params = tstnep_params;

end