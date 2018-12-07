function problem = tepsis_test_array(array_id)
% This is the HPC experiment for comparing different solution methods for
% the Transmission Explansion Problem with a Single Investment Stage TEPSIS
% The test array attempts to optimize the TEP with three sampling methods:
% Random Samples, K-means, and Latent Factor Approximation, and two
% optimization methods: MIP, and Benders Decomposition. After optimizing,
% the array checks the accuracy of the found plans.


%History            
%Version    Date        Who     Summary
%1          06/07/2018  JesseB  Adapted from LF_array_call for one stage TEP test


%% Initialize Problem Data
problem.bus_n = 312;
problem.gen_n = 980;
problem.non_dispatch_gen = 332:980;
problem.load_growth = 1.3;
problem.existing_plan = 1:654;


%% Determine Array ID Task
case_list = [2 4 6];
case_id= case_list(array_id);

switch case_id
    case {1, 3, 5}
        % Mixed integer program optimization model
        problem.opt_mod = "MIP";
    case {2, 4, 6}
        % Bender's decomposition optimization model
        problem.opt_mod = "bender";
end

switch case_id
    case 1
        % latent factor approximation samples with MIP opt model
        problem.samp_method = "LF";
        samp_n = 75; 
        lf_data = matfile('lf_approx_data.mat');
        
        scen_list = lf_data.scen_list;
        scen_w = lf_data.scen_w_pos;
        psub_m = lf_data.scen_mean;
        global_mean = lf_data.global_mean;  
        
    case 2
       % latent factor approximation samples with benders decomposition opt model
        problem.samp_method = "LF";
        samp_n = 75;        
        lf_data = matfile('lf_approx_data.mat');
        
        scen_list = lf_data.scen_list;
        scen_w = lf_data.scen_w;
        psub_m = lf_data.scen_mean;
        global_mean = lf_data.global_mean;
     
    case {3, 4}
        % monte carlo samples
        problem.samp_method = "MC";
        samp_n = 25;
        mc_samps = matfile('mc_scenarios.mat');
        
        scen_list = mc_samps.scen_list;
        scen_w = mc_samps.scen_w;
        
    case {5, 6}
        % k means samples
        problem.samp_method = "KM";
        samp_n = 60;
        k_data = matfile('k_means_data.mat');
        
        scen_w = k_data.scen_w;
        scen_load = k_data.scen_load;
        scen_renew = k_data.scen_renew;
        scen_varcost = k_data.scen_varcost;
end



%% Storage Initialization
outfile_name = sprintf('%s_%d','array',array_id);
output = matfile(outfile_name, 'Writable',true);
output.obj_val = zeros(samp_n,1);
output.lower_bound = zeros(samp_n,1);
output.opt_gap = zeros(samp_n,1);
output.run_time = zeros(samp_n,1);

lines_built = cell(samp_n,1);

%% Run Optimizations
for s_idx = 1:samp_n
    % write problem structure 
    switch case_id
        case {1,2}
        % set latent factor problem data    
        problem.scen_list = scen_list{s_idx}';
        problem.scen_w = scen_w{s_idx};
        problem.psub_m = psub_m{s_idx};
        problem.global_mean = global_mean(s_idx,1);         
        
        case {3,4}
        % set monte carlo problem data    
        problem.scen_list = scen_list{s_idx};
        problem.scen_w = scen_w{s_idx};    
        
        case {5,6}
        % set k means problem data
        problem.scen_load = scen_load{s_idx};
        problem.scen_renew = scen_renew{s_idx};
        problem.scen_VarCost = scen_varcost{s_idx};
        problem.scen_w = scen_w{s_idx};    
        problem.scen_list = 1:length(scen_w{s_idx});          
    
    end
    
    % initialize values for error checking
    problem.tot_cost = -1; problem.lower_bound = -1; problem.opt_gap = -1; problem.run_time = -1;
    problem.lines = 0;
    
    % run optimization
    problem = run_tep_test_gams(problem);
    
    % save output for run
    output.obj_val(s_idx,1)= problem.tot_cost;
    output.lower_bound(s_idx,1)= problem.lower_bound;
    output.opt_gap(s_idx,1)= problem.opt_gap;
    output.run_time(s_idx,1)= problem.run_time;
    lines_built{s_idx}= problem.lines;
    output.lines_built = lines_built; 
end

output.lines_built = lines_built;

end