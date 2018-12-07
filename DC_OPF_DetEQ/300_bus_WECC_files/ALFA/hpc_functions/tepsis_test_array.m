function problem = tepsis_test_array(array_id, case_id)
% This is the HPC experiment for comparing different solution methods for
% the Transmission Explansion Problem with a Single Investment Stage TEPSIS
% The test array attempts to optimize the TEP with three sampling methods:
% Random Samples, K-means, and Latent Factor Approximation, and two
% optimization methods: MIP, and Benders Decomposition. After optimizing,
% the array checks the accuracy of the found plans.


%History            
%Version    Date        Who     Summary
%1          06/07/2018  JesseB  Adapted from LF_array_call for one stage TEP test
%2          08/10/2018  JesseB  Added Importance Sampling
%3          08/13/2018  JesseB  Reworked for second experiment


%% Initialize Problem Data
problem.bus_n = 312;
problem.gen_n = 980;
problem.non_dispatch_gen = 332:980;
problem.load_growth = 1.3;
problem.existing_plan = 1:654;
run_time_list = 3600*[.25 .5 .75 1 1.5 2 3 4];

%% Determine Array ID Task
%case_list = [2 4 6];
%case_id= case_list(array_id);
%case_id = 5;
switch case_id
    case {1, 2, 3, 5, 7}
        % Mixed integer program optimization model
        problem.opt_mod = "MIP";
    case {4, 6}
        % Bender's decomposition optimization model
        problem.opt_mod = "bender";
end




switch case_id
    case 1
        % latent factor approximation samples with MIP opt model
        problem.samp_method = "LF";
        outfile_string = "LFS";
        lf_data = matfile('lf_10_data.mat');
        
        samp_n = length(lf_data.global_mean);
        scen_list = lf_data.scen_list;
        scen_w = lf_data.scen_w_pos;
        psub_m = lf_data.scen_mean;
        global_mean = lf_data.global_mean;  
        
    case 2
       % latent factor approximation samples with benders decomposition opt model
        problem.samp_method = "LF";
        outfile_string = "LFB";   
        lf_data = matfile('lf_30_data.mat');
        
        samp_n = length(lf_data.global_mean);
        scen_list = lf_data.scen_list;
        scen_w = lf_data.scen_w_pos;
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
        outfile_string = "KM";
        k_data = matfile('km_data.mat');
        samp_n = length(k_data.cluster_list);
        
        scen_n = k_data.cluster_list;
        scen_count = k_data.scen_w;
        scen_load = k_data.scen_load;
        scen_renew = k_data.scen_renew;
        scen_varcost = k_data.scen_varcost;
              
    case 7
        problem.samp_method = "IS";
        outfile_string = "IS";
        is_data = matfile('is_data.mat');
        
        samp_n = length(is_data.scen_list);
        scen_list = is_data.scen_list;
        scen_w = is_data.scen_w;     
end

samp_n = 8;
base_idx = (array_id-1) + 40;

%% Storage Initialization
outfile_name = sprintf('%s_%d',outfile_string, array_id);
output = matfile(outfile_name, 'Writable',true);
output.obj_val = zeros(samp_n,1);
output.lower_bound = zeros(samp_n,1);
output.opt_gap = zeros(samp_n,1);
output.run_time = zeros(samp_n,1);
output.run_time_limit = run_time_list;


lines_built = cell(samp_n,1);

%% Run Optimizations
for s_idx = 1:samp_n
    % write problem structure 
    switch case_id
        case {1,2}
        % set latent factor problem data    
        problem.scen_list = scen_list{base_idx}';
        problem.scen_w = scen_w{base_idx};
        problem.psub_m = psub_m{base_idx};
        problem.global_mean = global_mean(base_idx);         
        
        
        case {3,4}
        % set monte carlo problem data    
        problem.scen_list = scen_list{base_idx + s_idx};
        problem.scen_w = scen_w{base_idx + s_idx};    
        
        case {5,6}
        % set k means problem data
        problem.scen_load = scen_load{base_idx + s_idx};
        problem.scen_renew = scen_renew{base_idx + s_idx};
        problem.scen_VarCost = scen_varcost{base_idx + s_idx};
        problem.scen_w = scen_count{base_idx + s_idx}./8736;    
        problem.scen_list = 1:scen_n(base_idx + s_idx); 
        
        case 7
        problem.scen_list = scen_list{base_idx};
        problem.scen_w = scen_w{base_idx};

    end
    
    % initialize values for error checking
    problem.tot_cost = -1; problem.lower_bound = -1; problem.opt_gap = -1; problem.run_time = -1;
    problem.lines = 0;
    
    % run optimization
    problem.run_time_lim = run_time_list(s_idx);
    problem = run_tep_test_gams(problem);
    
    % save output for run
    output.obj_val(s_idx,1)= problem.tot_cost;
    output.lower_bound(s_idx,1)= problem.lower_bound;
    output.opt_gap(s_idx,1)= problem.opt_gap;
    output.run_time(s_idx,1)= problem.run_time;
    lines_built{s_idx}= problem.lines;
    output.lines_built = lines_built; 
end

%output.lines_built = lines_built;

end