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

%% Organize Candidate Plan for Optimal Operation Cost Run
cand_n = sum(params.line.cand);
built_n = sum(params.line.built);
dec_built = logical(ones(built_n, cand_n));
dec_new = logical(eye(cand_n));
dec_lines = [dec_built; dec_new];
dec_lines = [logical(params.line.y), dec_lines];
cand_cost = cell(cand_n,1);
problem.init_time = toc;
tic;
for c_idx = 1:cand_n
params.candidate.imp = params.line.imp(dec_lines(:,c_idx));
params.candidate.res = params.line.res(dec_lines(:,c_idx));
params.candidate.from_to = params.line.from_to(dec_lines(:,c_idx),:);
params.candidate.max_flow = params.line.max_flow(dec_lines(:,c_idx));
params.candidate.n = sum(params.line.y);
params.candidate.loss.const = params.line.loss.const;
params.candidate.loss.slope = params.line.loss.slope;
params.candidate.loss.n = size(params.candidate.loss.slope,1);


%% Run Solution in Parallel
op_cost = cell(params.scen.n,1);
for scen = 1:params.scen.n
    op_cost{scen} = DC_OPF_operations(params, dec_lines(:,c_idx), scen);
end

%% Solution Output
op_cost = cell2mat(op_cost);
cand_cost{c_idx} = params.scen.p'*op_cost;

end
problem.runtime = toc;
problem.cand_cost = cell2mat(cand_cost);
problem.solution_value = min(problem.cand_cost);
problem.params = params;

end