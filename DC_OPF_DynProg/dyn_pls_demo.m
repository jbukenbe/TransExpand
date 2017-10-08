function problem = dyn_pls_demo()
% This function demonstrates the use of pls regression to sample candidate
% plans for the transmission expansion problem and make and informed search
% for new plans

%History            
%Version    Date        Who     Summary
%1          10/07/2017  JesseB  Adapted from dyn_DC_OPF_demo
%2          10/07/2017  JesseB  Implimented initial sample and pls search


% To Do: separate optimization parameters from pls parameters      

%% Initialize Data
tic;
params = DC_OPF_init;
y_var = [51 52 53 55 57 60 62 64 68 69 72 75 82 92 94 96 98 100 118 131]';
%y_var = [51 52 53 55 57 72 75 82 118 131]';

%% Organize Candidate Plan for Optimal Operation Cost Run
%calculate number of candidate lines and plans
cand_n = sum(y_var>0);
params.cand.n = cand_n;
plan_n = 2^cand_n;
built_n = sum(params.line.built);
dec_built = logical(params.line.built);
scen_op_cost = zeros(plan_n, params.scen.n);
cand_cost = cell(plan_n,1);
plan_id = cell(plan_n,1);
problem.init_time = toc;
tic;

%make initial sample of plans
samp_id = randperm(plan_n,params.initial_samp_n)';
cur_best_val = 10000;
err = 1;
while err > params.stop_err
    for c_idx = 1:size(samp_id,1)
        new_line_idx = y_var.*de2bi(samp_id(c_idx)-1,cand_n)';
        new_line_idx = nonzeros(new_line_idx);
        dec_lines = dec_built;
        dec_lines(new_line_idx) = true;
        params.cand.imp = params.line.imp(dec_lines);
        params.cand.res = params.line.res(dec_lines);
        params.cand.from_to = params.line.from_to(dec_lines,:);
        params.cand.max_flow = params.line.max_flow(dec_lines);
        
        
        %% Run Solution in Parallel
        op_cost = cell(params.scen.n,1);
        parfor scen = 1:params.scen.n
            op_cost{scen} = DC_OPF_operations(params, dec_lines, scen);
        end
        
        %% Solution Output
        scen_op_cost(samp_id(c_idx),:) = cell2mat(op_cost)';
        cand_cost{samp_id(c_idx)} = scen_op_cost(samp_id(c_idx),:)*params.scen.p;
        plan_id{samp_id(c_idx)} = samp_id(c_idx);
    end
problem.plan_id = cell2mat(plan_id);
problem.scen_op_cost = scen_op_cost;
problem.cand_cost = cell2mat(cand_cost);
[new_best_val,best_plan_id] = min(problem.cand_cost);
err = (cur_best_val - new_best_val)/new_best_val;
cur_best_val = new_best_val;
if err > params.stop_err
    problem.params = params;
    problem = pls_val_est(problem);
    samp_id = problem.new_plan_id;
end
end

problem.solution_value = new_best_val;
problem.solution = nonzeros(de2bi(problem.plan_id(best_plan_id)-1,cand_n)'.*y_var);
problem.params = params;
problem.runtime = toc;
end