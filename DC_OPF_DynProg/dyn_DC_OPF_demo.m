function problem = dyn_DC_OPF_demo()
% This shell program runs multiple scenarios for the DC OPF problem for a
% number of candidate plans (Currently only one candidate) and weights the
% solution based on the scenario probability of occurence

%History            
%Version    Date        Who     Summary
%1          09/10/2017  JesseB  Initial Version
%2          11/02/2017  JesseB  Remodeled from dyn_pls_demo for parallelization

%% Initialize Data
% start initialization timer
tic;
% read input data files
params = DC_OPF_init;

% select candidate lines to include in analysis TODO
%y_var = [51 52 53 55 57 60 62 63 64 65 68 69 72 75 82 84 85 86 92 94 96 98 100 118 131]';
%y_var = [51 52 53 55 57 60 62 64 68 69 72 75 82 92 94 96 98 100 118 131]';
%y_var = [51 52 53 55 57 60 62 64 69 72 75 82 94 118 131]';
y_var = [51 52 53 55 57 72 75 82 118 131]';

% load relavant line costs
params.new_line_cost = params.line.cost(y_var);


%% Organize Candidate Plan for Optimal Operation Cost Run
% calculate number of candidate lines and plans
cand_n = sum(y_var>0);
params.cand.n = cand_n;
plan_n = 2^cand_n;
dec_built = logical(params.line.built);

% create parallel cell variables
par_scen_op_cost = cell(plan_n, params.scen.n);
par_cand_op_cost = cell(plan_n, 1);
par_cand_full_cost = cell(plan_n,1);

problem.init_time = toc;
tic;
    
% look at all lines
parfor c_idx = 1:plan_n
% load relavent data for DC OPF for this plan
    new_line_idx = y_var.*de2bi(c_idx-1,cand_n)';
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
end

problem.scen_op_cost = cell2mat(par_scen_op_cost); clear par_scen_op_cost;
problem.cand_op_cost = cell2mat(par_cand_op_cost); clear par_cand_op_cost;
problem.cand_full_cost = cell2mat(par_cand_full_cost); clear par_full_cost;
problem.plan_id = [1:plan_n]';

%% Output Data
[problem.solution_value, problem.solution] = min(problem.cand_full_cost);
problem.solution = problem.solution-1;
problem.params = params;
problem.runtime = toc;
write_struct(problem);

end