function status = LF_launchpad()
% This is the launchpad program to test the latent factor stocastic
% approximator algorithm on HPC

%History            
%Version    Date        Who     Summary
%1          02/19/2018  JesseB  Initial Version

% To Do:
%   Allow for easy switching between models and model sizes
%   Store data for Crossvalidation set
%   Develop Crossvalidation sets for different problems
%   Finish error quantification metrics
%   Sampling method for initial set



%% Initializaion
status = 0;
test_rep_n = 2;
initial_samp_list = 20:20:200;
lat_fact_list = 1:4;
samp_per_lf_list = 1;

init_samp_n = length(initial_samp_list);
lat_fact_n = length(lat_fact_list);
samp_per_lf_n = length(samp_per_lf_list);

% Run crossvalidation set
cross_val_prob.samp_n = 100; 
cross_val_prob = LF_crossval_set(cross_val_prob, 'calculate');

% Create Data log matricies
%initial_samp_runtime_log = zeros(test_rep_n, init_samp_n, lat_fact_n, samp_per_lf_n);
%second_samp_runtime_log = zeros(test_rep_n, init_samp_n, lat_fact_n, samp_per_lf_n);

%observed_cost_log = zeros(test_rep_n, init_samp_n, lat_fact_n, samp_per_lf_n);
%approx_cost_log = zeros(test_rep_n, init_samp_n, lat_fact_n, samp_per_lf_n);
err_ss_log = zeros(test_rep_n, init_samp_n, lat_fact_n, samp_per_lf_n);
approx_r_sq_log = zeros(test_rep_n, init_samp_n, lat_fact_n, samp_per_lf_n);
observed_r_sq_log = zeros(test_rep_n, init_samp_n, lat_fact_n, samp_per_lf_n);


%% Run Iterations
for n_idx = 1:test_rep_n
    for i_idx = 1:init_samp_n
% Make initial sample for speed and comparability
        init_samp = initial_samp_list(i_idx);    
        test_prob = LF_test_init(init_samp);

        for f_idx = 1:lat_fact_n
% Number of Latent Factors
            test_prob.LF_n = lat_fact_list(f_idx);

            for s_idx = 1:samp_per_lf_n
% Scenarios per Latent Factor
                test_prob.samp_per = samp_per_lf_list(s_idx);            

% Run Tests
                test_prob = svd_approx(test_prob, cross_val_prob);           

% Record data
                %initial_samp_runtime_log(n_idx, i_idx, f_idx, s_idx) = test_prob.init_samp_runtime;
                %second_samp_runtime_log(n_idx, i_idx, f_idx, s_idx) = test_prob.second_samp_runtime;

                %approx_cost_log(n_idx, i_idx, f_idx, s_idx) = test_prob.approx_cost;
                %observed_cost_log(n_idx, i_idx, f_idx, s_idx) = test_prob.actual_cost;
                err_ss_log(n_idx, i_idx, f_idx, s_idx) = test_prob.err_sum_squares;
                approx_r_sq_log(n_idx, i_idx, f_idx, s_idx) = test_prob.r_squared_est;
                observed_r_sq_log(n_idx, i_idx, f_idx, s_idx) = test_prob.obs_r_squared;
            end
        end
    end
end

%% Write recorded data to output file
%log.initial_samp_runtime = initial_samp_runtime_log;
%log.second_samp_runtime = second_samp_runtime_log;
%log.total_runtime = total_runtime_log;
%log.approx_cost = approx_cost_log;
%log.observed_cost = observed_cost_log;
log.err_ss = err_ss_log;
log.approx_r_sq = approx_r_sq_log;
log.observed_r_sq = observed_r_sq_log;

filename = 'LF_logs.mat';
output = matfile(filename);
output.log = log;
status = 1
end
