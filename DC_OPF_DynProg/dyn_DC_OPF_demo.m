function problem = dyn_DC_OPF_demo()
% This shell program runs multiple scenarios for the DC OPF problem for a
% number of candidate plans (Currently only one candidate) and weights the
% solution based on the scenario probability of occurence

%History            
%Version    Date        Who     Summary
%1          09/10/2017  JesseB  Initial Version

%% Initialize Data
tic;
params = DC_OPF_init;
y_var = [51 52 53 55 57 72 75 82 118 131]';

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
for c_idx = 1:plan_n
new_line_idx = y_var.*de2bi(c_idx-1,cand_n)';
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
scen_op_cost(c_idx,:) = cell2mat(op_cost)';
cand_cost{c_idx} = scen_op_cost(c_idx,:)*params.scen.p;
plan_id{c_idx} = c_idx;

end
problem.runtime = toc;
problem.plan_id = cell2mat(plan_id);
problem.scen_op_cost = scen_op_cost;
problem.cand_cost = cell2mat(cand_cost);
problem.solution_value = min(problem.cand_cost);
problem.params = params;

end