function problem = dyn_pls_demo()
% This function demonstrates the use of pls regression to sample candidate
% plans for the transmission expansion problem and make an informed search
% for new plans

%History            
%Version    Date        Who     Summary
%1          10/07/2017  JesseB  Adapted from dyn_DC_OPF_demo
%2          10/07/2017  JesseB  Implimented initial sample and pls search
%3          10/08/2017  JesseB  Added line costs
%4          10/08/2017  JesseB  Updated for parallel plan solutions
%5          10/09/2017  JesseB  Fractional Factorial initial sample
%6          11/15/2017  JesseB  Removed parallel operations to par_plan_ops
%7          11/15/2017  JesseB  Streamlined for 2 phase search algorithm
%8          12/02/2017  JesseB  Moved candidate line data to initialization


%% Initialize Data
% set problem run list
problem_size_run_list = [196, 50, 45,40,35,30];

% read input data files
for run_idx = 1:1
init_time = tic;
problem.z_idx = 1;  z_idx=1;
problem.problem_size = problem_size_run_list(run_idx);
problem.params = DC_OPF_init(problem);    
problem.lock_on{1} = [];
problem.lock_off{1} = [];

%% Organize Candidate Plan for Optimal Operation Cost Run
% calculate number of candidate lines and plans
line_id = problem.params.cand.line_id{1};
cand_n = problem.params.cand.n;

% make initial sample of plans
%problem.samp_id = randperm(plan_n,params.initial_samp_n)';
problem.samp = fraction_fact_samp(problem);
problem.plan{1} = problem.samp;
problem.samp_range(1,problem.z_idx) = 1;
problem.samp_range(2,problem.z_idx) = size(problem.samp,1); 
problem.init_time = toc(init_time);
run_time = tic;

%% Run DC OPF for Initial Sample
problem = par_plan_ops(problem);
problem.run_time = toc(run_time);
[best_plan_full_cost, best_plan_id] = min(problem.cand_full_cost);
best_plan = problem.params.cand.line_id{z_idx}.*double(nonzeros(problem.plan{z_idx}(best_plan_id,:)))';
best_set = 1;
fprintf('%d%s%6.2f\n',problem.z_idx,'  Current best solution: ',best_plan_full_cost);

%% Save SVD Latent Factor Info
problem.params.svd.use_latent_fac = 1;
[~,problem.params.svd.s_values, problem.params.svd.directions] = svd(problem.scen_op_cost);
[~,factor_id] =sort(abs(problem.params.svd.directions(:,1:problem.params.svd.scen_n)),'descend');
factor_id = unique(factor_id','rows','stable'); 
problem.params.svd.latent_scen = factor_id(1:problem.params.svd.scen_n);
problem.params.svd.filler_scen = setdiff(1:problem.params.scen.n,problem.params.svd.latent_scen);
problem.params.svd.good_plan_threshold = min(problem.cand_full_cost);

%% Loop Through PLS fits
while problem.z_idx < 10  
    %% Update SVD info
    problem.params.svd.good_plan_threshold = min(problem.cand_op_cost)*1.05;
    
    %% Run PLS Regression
    reg_time = tic;
    problem = pls_val_est(problem);
    problem.pls_reg_time(z_idx) = toc(reg_time);
    find_samp = tic;
    
    problem.z_idx = problem.z_idx+1;    z_idx = problem.z_idx;
    
    %% Get Refined Sample
    problem = new_samp_from_beta(problem);
    problem.plan{problem.z_idx} = problem.samp;
    problem.samp_range(1,problem.z_idx) = problem.samp_range(2,problem.z_idx-1)+1;
    problem.samp_range(2,problem.z_idx) = problem.samp_range(1,problem.z_idx)+ size(problem.samp,1)-1; 
    problem.find_refined_samp_time(z_idx-1) = toc(find_samp);
    run_time = tic;

    %% Run DC OPF for Refined Sample
    % Save data from first sample
    saved_scen_op_cost = problem.scen_op_cost; clear problem.scen_op_cost;
    saved_cand_op_cost = problem.cand_op_cost; clear problem.cand_op_cost;
    saved_cand_full_cost = problem.cand_full_cost; clear problem.cand_full_cost;

    % Run refined sample
    problem = par_plan_ops(problem);
    [best_set_plan_full_cost, best_set_plan] = min(problem.cand_full_cost);
    itr_best_plan = nonzeros(problem.params.cand.line_id{z_idx}.*double((problem.plan{z_idx}(best_set_plan,:)))');
    itr_best_plan = [itr_best_plan;problem.lock_on{z_idx};];
    
    % Reload and combine sample data
    problem.scen_op_cost = [saved_scen_op_cost; problem.scen_op_cost]; clear saved_scen_op_cost;
    problem.cand_op_cost = [saved_cand_op_cost; problem.cand_op_cost]; clear saved_cand_op_cost;
    problem.cand_full_cost= [saved_cand_full_cost; problem.cand_full_cost]; clear saved_cand_full_cost;

    % Find the best plan and its costs
    [~, best_plan_id] = min(problem.cand_full_cost);
    if best_set_plan_full_cost < best_plan_full_cost
        best_plan_full_cost = best_set_plan_full_cost;
        best_set = z_idx;
        best_plan = itr_best_plan;
        fprintf('%d%s%6.2f\n',problem.z_idx,'  New solution found: ',best_set_plan_full_cost);
    else
        fprintf('%d%s%6.2f\n',problem.z_idx,'  This best solution: ',best_set_plan_full_cost);
    end
    problem.run_time(z_idx) = toc(run_time);
end

%% Output Data
problem.solution_value = best_plan_full_cost;
problem.solution_id = best_plan_id;
problem.best_set = best_set;
problem.solution_lines = best_plan;
problem.runtime = toc(init_time);
filename = sprintf('%s_%d','output',run_idx);
write_struct(problem, filename);
clear problem
end
end