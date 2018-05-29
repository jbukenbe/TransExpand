
samp_size_list = [10 30 50];
lf_size_list = [2 3 4 5 6 8 10 15 20 30 40 50 60 70 80 90 100];

runs = 20;
samp_length = length(samp_size_list);
lf_length = length(lf_size_list);
log_ape = cell(samp_length,lf_length,runs);
log_ess = cell(samp_length,lf_length,runs);
log_tot = cell(samp_length,lf_length,runs);
log_rse = cell(samp_length,lf_length,runs);
log_rso = cell(samp_length,lf_length,runs);
log_prst = cell(samp_length,lf_length,runs);
log_mc_tot = cell(samp_length,lf_length,runs);
for t_idx = 1:runs
prob.samp_per = 1;
prob.params.scen.n = 8736;
subset = [scen_op_cost(100:163,:);output];
cross_val_prob.scen_op_cost = output;
prob.scen_op_cost = [];
for samp_idx = 1:samp_length
    samp_size = samp_size_list(samp_idx);
    if samp_idx == 1
        prob.scen_op_cost = subset(randperm(364,samp_size),:);
        last_samp_size = samp_size;
    else 
        prob.scen_op_cost = [prob.scen_op_cost;subset(randperm(364,samp_size-last_samp_size),:)];
    end
%% Initialize approximation tools for future use if needed
% mean center scenario data
prob.mean_scen_cost = mean(prob.scen_op_cost);
A = prob.scen_op_cost - prob.mean_scen_cost;
% select scenarios from svd latent factors
[U, S, V] = svd(A);
prob.S = S;
prob.V = V;
    
    for lf_idx = 1:lf_length
        if lf_size_list(lf_idx) <= samp_size
            par_prob = prob;
            par_prob.LF_n = lf_size_list(lf_idx);
        else
            par_prob = prob;
            par_prob.LF_n = samp_size;
            par_prob.samp_per = lf_size_list(lf_idx)/samp_size;
        end           
            
        par_prob = svd_approx(par_prob, cross_val_prob);
        log_ape{samp_idx, lf_idx,t_idx} = par_prob.mean_abs_percent_err; 
        log_ess{samp_idx, lf_idx,t_idx} = par_prob.err_sum_squares;
        log_tot{samp_idx, lf_idx,t_idx} = par_prob.tot_abs_percent_err;
        log_rse{samp_idx, lf_idx,t_idx} = par_prob.r_squared_est;
        log_rso{samp_idx, lf_idx,t_idx} = par_prob.obs_r_sq2;
        log_prst{samp_idx, lf_idx,t_idx} = par_prob.tot_r_sq2;
        log_mc_tot{samp_idx, lf_idx,t_idx} = par_prob.mc_tot_percent_err;
    end
end
end
log_a = cell2mat(log_ape);
log_e = cell2mat(log_ess);
log_ro = cell2mat(log_rso);
log_re = cell2mat(log_rse);
log_tr = cell2mat(log_prst);
log_t = cell2mat(log_tot);
log_tmc = cell2mat(log_mc_tot);
%% Load Data

%{
load('20_line_output.mat');
problem.params = DC_OPF_init;

y_var = [51 52 53 55 57 60 62 64 68 69 72 75 82 92 94 96 98 100 118 131]';
%y_var = [51 52 53 55 57 60 62 64 69 72 75 82 94 118 131]';
%y_var = [51 52 53 55 57 72 75 82 118 131]';

%% Initialization
% PLS regression initialization
problem.params.new_line_cost = problem.params.line.cost(y_var);
problem.cand_full_cost = cand_full_cost;
%problem.real_op = cand_op_cost;
problem.params.interaction = 1;
problem.params.cand.n = length(y_var);

% full problem initialization
full_n = size(plan_id,1);
best_val = min(cand_full_cost);
[~,real_rank] = sort(cand_full_cost);
real_top_set = real_rank(1:100);

% sample and run init
samp_n_list = [100:50:500]';
samp_range_n = size(samp_n_list,1);
k_runs = 100;

% Output Storage init
avg_found_best = zeros(1,samp_range_n);
avg_best_rank= zeros(1,samp_range_n);
avg_best_qual= zeros(1,samp_range_n);
avg_prop_best= zeros(1,samp_range_n);


%% Run and Save Results
% iterate through sample size
for samp_idx = 1:samp_range_n
% set sample size
    init_samp_n = samp_n_list(samp_idx);
    
% run storage init
    found_best = cell(k_runs,1);
    best_rank = cell(k_runs,1);
    best_qual = cell(k_runs,1);
    prop_best = cell(k_runs,1);
    
% run k searches
    for run_idx = 1:k_runs
% make sample
        samp_plan_id = randperm(full_n, init_samp_n)';
        plan_op_val = cand_op_cost(samp_plan_id);

% make paralell problem structure
        par_prob = problem;
        par_prob.plan_id = samp_plan_id;
        par_prob.cand_op_cost = plan_op_val;
        
% run search
        [best_found, best_found_id] = pls_val_est(par_prob);
        
% store search data for run
        found_best{run_idx} = (best_found(1) == best_val);
        best_rank{run_idx} = find(real_rank == best_found_id(1),1);
        best_qual{run_idx} = best_found(1);
        prop_best{run_idx} = (100-length(setdiff(best_found_id,real_top_set)));    
    end
    
% save average run data in output files
    avg_found_best(samp_idx) = mean(cell2mat(found_best)); clear found_best;
    avg_best_rank(samp_idx) = mean(cell2mat(best_rank)); clear best_rank;
    avg_best_qual(samp_idx) = mean(cell2mat(best_qual)); clear best_qual;
    avg_prop_best(samp_idx) = mean(cell2mat(prop_best)); clear prop_best;
end
figure 
subplot (2,2,1)
plot(samp_n_list, avg_found_best);
title('Frequency of Finding Best Plan'); xlabel('Initial Sample n'); ylabel('Frequency');

subplot (2,2,2)
plot(samp_n_list, avg_prop_best);
title('Proportion of 100 Best Plans Found'); xlabel('Initial Sample n'); ylabel('%');

subplot (2,2,3)
plot(samp_n_list, avg_best_rank);
title('Avg Rank of Best Found Plan'); xlabel('Initial Sample n'); ylabel('Rank');

subplot (2,2,4)
plot(samp_n_list, avg_best_qual);
title('Avg Cost of Best Found Plan'); xlabel('Initial Sample n'); ylabel('Cost');
%}

