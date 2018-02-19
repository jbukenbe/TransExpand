function status = LF_launchpad()
% This is the launchpad program to test the latent factor stocastic
% approximator algorithm on HPC

%History            
%Version    Date        Who     Summary
%1          02/19/2018  JesseB  Initial Version


%% Initializaion
status = 0;
test_rep_n = 2;
initial_samp_list = 20:20:200;
lat_fact_list = 1:20;
samp_per_lf = 1:3;

init_samp_n = length(initial_samp_list);
lat_fact_n = length(lat_fact_list);
samp_per_lf_n = length(samp_per_lf);

initial_samp_runtime_log = zeros(test_rep_n, init_samp_n, lat_fact_n, samp_per_lf_n);
second_samp_runtime_log = zeros(test_rep_n, init_samp_n, lat_fact_n, samp_per_lf_n);
total_runtime_log = zeros(test_rep_n, init_samp_n, lat_fact_n, samp_per_lf_n);

observed_cost_log = zeros(test_rep_n, init_samp_n, lat_fact_n, samp_per_lf_n);
approx_cost_log = zeros(test_rep_n, init_samp_n, lat_fact_n, samp_per_lf_n);
approx_ss_log = zeros(test_rep_n, init_samp_n, lat_fact_n, samp_per_lf_n);
observed_ss_log = zeros(test_rep_n, init_samp_n, lat_fact_n, samp_per_lf_n);
approx_r_sq_log = zeros(test_rep_n, init_samp_n, lat_fact_n, samp_per_lf_n);
observed_r_sq_log = zeros(test_rep_n, init_samp_n, lat_fact_n, samp_per_lf_n);


%% Run Iterations
for i_idx = 1:init_samp_n
% Initial Sample Size
    init_samp = init_samp_list(i_idx);    
    
    for f_idx = 1:lat_fact_n
% Number of Latent Factors
        lat_fact = lat_fact_list(f_idx);
        
        for s_idx = 1:samp_per_lf_n
% Scenarios per Latent Factor
            samp_per = samp_per_lf_list(s_idx);            
            
            for n_idx = 1:test_rep_n
% Run Tests
                test_prob = LF_test_init(init_samp, lat_fact, samp_per);
                test_prob = svd_approx(test_prob, 'initialize');
                test_prob = LF_test_secondary(test_prob);
                test_prob = svd_approx(test_prob, 'estimate');

% Record data
                initial_samp_runtime_log(n_idx, i_idx, f_idx, s_idx) = test_prob.init_samp_runtime;
                second_samp_runtime_log(n_idx, i_idx, f_idx, s_idx) = test_prob.second_samp_runtime;
                total_runtime_log(n_idx, i_idx, f_idx, s_idx) = test_prob.total_runtime;

                approx_cost_log(n_idx, i_idx, f_idx, s_idx) = test_prob.approx_cost;
                observed_cost_log(n_idx, i_idx, f_idx, s_idx) = test_prob.actual_cost;
                approx_ss_log(n_idx, i_idx, f_idx, s_idx) = test_prob.sum_squares_est;
                observed_ss_log(n_idx, i_idx, f_idx, s_idx) = test_prob.obs_sum_squares;
                approx_r_sq_log(n_idx, i_idx, f_idx, s_idx) = test_prob.r_squared_est;
                observed_r_sq_log(n_idx, i_idx, f_idx, s_idx) = test_prob.obs_r_squaree;
            end
        end
    end
end

%% Write recorded data to output file
log.initial_samp_runtime = initial_samp_runtime_log;
log.second_samp_runtime = second_samp_runtime_log;
log.total_runtime = total_runtime_log;
log.approx_cost = approx_cost_log;
log.observed_cost = observed_cost_log;
log.approx_ss = approx_ss_log;
log.observed_ss = observed_ss_log;
log.approx_r_sq = approx_r_sq_log;
log.observed_r_sq = observed_r_sq_log;

filename = 'LF_logs.mat';
output = matfile(filename);
output.log = log;
status = 1
end



%{
robust_log = zeros(20,20);
expected_cost_err_log = zeros(20,20);

for svd_scen = 1:20
    for init_idx = 10:10:200
        test_prob = dyn_pls_demo(saved_problem, init_idx+25,svd_scen);
        second_samp_approx = test_prob.scen_op_cost((26+init_idx):end,:);
        
        robust_error = (abs(max(saved_problem.second_samp_real(36:end,:)')-max(second_samp_approx')))./max(saved_problem.second_samp_real(36:end,:)');
        robust_mse = sum(robust_error)/3;
        robust_log(init_idx/10,svd_scen) = robust_mse;
                
        
        second_samp_err = (second_samp_approx - saved_problem.second_samp_real(36:end,:));
        second_samp_sq_err = second_samp_err.*second_samp_err;
        second_samp_mse = sum(sum(second_samp_sq_err))/(300*103);
        expected_cost_err_log(init_idx/10,svd_scen) = second_samp_mse;
    end
end

%}