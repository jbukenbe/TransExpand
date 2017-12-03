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

% To Do: separate optimization parameters from pls parameters      
%           Stop hardcoded line inclusion

%% Initialize Data
% set problem run list
problem_size_run_list = [50,45,40,35,30];

% read input data files
for run_idx = 1:5
init_time = tic;
problem_size = problem_size_run_list(run_idx);
params = DC_OPF_init(problem_size);    

%% Organize Candidate Plan for Optimal Operation Cost Run
% calculate number of candidate lines and plans
line_id = params.cand.line_id;
cand_n = params.cand.n;
plan_n = params.plan.n;

% make initial sample of plans
%problem.samp_id = randperm(plan_n,params.initial_samp_n)';
problem.samp_id = fraction_fact_samp(params);

% Move all data to problem structure
problem.params = params;
problem.init_time = toc(init_time);
problem.init_samp_run_time = tic;

%% Run DC OPF for Initial Sample
problem = par_plan_ops(problem);
problem.init_samp_run_time = toc(problem.init_samp_run_time);
problem.pls_reg_time = tic;

%% Run PLS Regression
problem = pls_val_est(problem);
problem.pls_reg_time = toc(problem.pls_reg_time);
problem.refined_samp_run_time = tic;

%% Run DC OPF for Refined Sample
% Save data from first sample
init_scen_op_cost = problem.scen_op_cost; clear problem.scen_op_cost;
init_cand_op_cost = problem.cand_op_cost; clear problem.cand_op_cost;
init_cand_full_cost = problem.cand_full_cost; clear problem.cand_full_cost;
problem.init_plan_id = problem.plan_id; clear problem.plan_id;

% Run refined sample 
problem = par_plan_ops(problem);

% Reload and combine sample data
problem.scen_op_cost = [init_scen_op_cost; problem.scen_op_cost]; clear init_scen_op_cost;
problem.cand_op_cost = [init_cand_op_cost; problem.cand_op_cost]; clear init_cand_op_cost;
problem.cand_full_cost= [init_cand_full_cost; problem.cand_full_cost]; clear init_cand_full_cost;
problem.refined_plan_id = problem.plan_id;
problem.plan_id = [problem.init_plan_id; problem.plan_id]; clear problem.plan_id

% Find the best plan and its costs
[best_plan_full_cost, best_plan_id] = min(problem.cand_full_cost);

problem.refined_samp_run_time = toc(problem.refined_samp_run_time);


%% Output Data
problem.solution_value = best_plan_full_cost;
problem.solution_id = best_plan_id;
problem.solution_lines = nonzeros(de2bi(problem.plan_id(best_plan_id)-1,cand_n)'.*line_id);
problem.runtime = toc(init_time);
filename = sprintf('%s_%d','output',run_idx);
write_struct(problem, filename);
clear problem params
end
end