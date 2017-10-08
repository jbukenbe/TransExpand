function problem = dyn_pls_demo()
% This function demonstrates the use of pls regression to sample candidate
% plans for the transmission expansion problem and make and informed search
% for new plans

%History            
%Version    Date        Who     Summary
%1          10/07/2017  JesseB  Adapted from dyn_DC_OPF_demo
%2          10/07/2017  JesseB  Implimented initial sample and pls search
%3          10/08/2017  JesseB  Added line costs
%4          10/08/2017  JesseB  Updated for parallel plan solutions

% To Do: separate optimization parameters from pls parameters      
%           Stop hardcoded line inclusion

%% Initialize Data
% start initialization timer
tic;
% read input data files
params = DC_OPF_init;

% select candidate lines to include in analysis TODO
y_var = [51 52 53 55 57 60 62 64 68 69 72 75 82 92 94 96 98 100 118 131]';
%y_var = [51 52 53 55 57 72 75 82 118 131]';

% load relavant line costs
params.new_line_cost = params.line.cost(y_var);


%% Organize Candidate Plan for Optimal Operation Cost Run
% calculate number of candidate lines and plans
cand_n = sum(y_var>0);
params.cand.n = cand_n;
plan_n = 2^cand_n;
dec_built = logical(params.line.built);
scen_op_cost = cell(plan_n, params.scen.n);
cand_op_cost = cell(plan_n, 1);
cand_full_cost = cell(plan_n,1);
plan_id = cell(plan_n,1);


% make initial sample of plans
samp_id = randperm(plan_n,params.initial_samp_n)';
samp_size = length(samp_id);
par_scen_op_cost = cell(samp_size, params.scen.n);
par_cand_op_cost = cell(samp_size, 1);
par_cand_full_cost = cell(samp_size,1);
par_plan_id = cell(samp_size,1);
cur_best_val = 10000;
err = 1;

problem.init_time = toc;
tic;
% begin search while improvement is less than threshold
while err > params.stop_err
    
% look at all sample lines
    parfor c_idx = 1:samp_size
% load relavent data for DC OPF for this plan
        new_line_idx = y_var.*de2bi(samp_id(c_idx)-1,cand_n)';
        new_line_idx = nonzeros(new_line_idx);
        dec_lines = dec_built;
        dec_lines(new_line_idx) = true;
        
% copy parameter data for parfor 
        par_params = params;
        par_params.cand.imp = params.line.imp(dec_lines);
        par_params.cand.res = params.line.res(dec_lines);
        par_params.cand.from_to = params.line.from_to(dec_lines,:);
        par_params.cand.max_flow = params.line.max_flow(dec_lines);
        
        
%% Run Scenario Optimization in Parallel
% run optimization to get operating cost of this plan for each scenario
        op_cost = cell(params.scen.n,1);
        for scen = 1:params.scen.n
            op_cost{scen} = DC_OPF_operations(par_params, dec_lines, scen);
        end
        
%% Parallel Data Cleanup
        par_scen_op_cost(c_idx,:) = op_cost';
        op_cost = cell2mat(op_cost);
        par_cand_op_cost{c_idx} = op_cost'*par_params.scen.p;
        par_cand_full_cost{c_idx} = par_cand_op_cost{c_idx} + sum(par_params.line.cost(new_line_idx));
        par_plan_id{c_idx} = samp_id(c_idx);
    end
    
scen_op_cost(samp_id,:) = par_scen_op_cost;
cand_op_cost(samp_id) = par_cand_op_cost;
cand_full_cost(samp_id) = par_cand_full_cost;
plan_id(samp_id) = par_plan_id;
    
    
%% Run PLS Search for new plans
% write data to problem structure to pass to pls estimator
    problem.scen_op_cost = cell2mat(scen_op_cost);
    problem.cand_op_cost = cell2mat(cand_op_cost);
    problem.cand_full_cost = cell2mat(cand_full_cost);
    problem.plan_id = cell2mat(plan_id);
    
% calculate improvement of plans
    [new_best_val, best_plan_id] = min(problem.cand_full_cost);
    err = (cur_best_val - new_best_val)/new_best_val;
    cur_best_val = new_best_val;
    
% run pls search if improvements are still being found
    if err > params.stop_err
        problem.params = params;
        problem = pls_val_est(problem);
        samp_id = problem.new_plan_id;
        samp_size = length(samp_id);
        par_scen_op_cost = cell(samp_size, params.scen.n);
        par_cand_op_cost = cell(samp_size, 1);
        par_cand_full_cost = cell(samp_size,1);
        par_plan_id = cell(samp_size,1);
    end
end

%% Output Data
problem.solution_value = new_best_val;
problem.solution = nonzeros(de2bi(problem.plan_id(best_plan_id)-1,cand_n)'.*y_var);
problem.params = params;
problem.runtime = toc;
end