%% Load Data
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


