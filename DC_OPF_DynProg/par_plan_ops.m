function problem = par_plan_ops(problem)
% This function takes the problem structure for the transmission expansion
% problem and solves the plan costs in parallel 

%History            
%Version    Date        Who     Summary
%1          10/07/2017  JesseB  Adapted from dyn_pls_demo

%% Load needed problem parameters
params = problem.params;
samp_id = problem.samp_id;
samp_size = length(samp_id);
cand_n = params.cand.n;
y_var = params.cand.line_id;
dec_built = params.line.dec_built;

%% Initialize Storage
par_scen_op_cost = cell(samp_size, params.scen.n);
par_cand_op_cost = cell(samp_size, 1);
par_cand_full_cost = cell(samp_size,1);
par_plan_id = cell(samp_size,1);

    
%% Run all sample plans in parallel
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
    op_cost = cell(par_params.scen.n,1);
    for scen = 1:par_params.scen.n
        op_cost{scen} = DC_OPF_operations(par_params, dec_lines, scen);
    end

%% Store Output Data
    par_scen_op_cost(c_idx,:) = op_cost';
    op_cost = cell2mat(op_cost);
    par_cand_op_cost{c_idx} = op_cost'*par_params.scen.p;
    par_cand_full_cost{c_idx} = par_cand_op_cost{c_idx} + sum(par_params.line.cost(new_line_idx));
    par_plan_id{c_idx} = samp_id(c_idx);
end
        
%% Parallel Data Cleanup
problem.scen_op_cost = cell2mat(par_scen_op_cost);
problem.cand_op_cost = cell2mat(par_cand_op_cost);
problem.cand_full_cost = cell2mat(par_cand_full_cost); 
problem.plan_id = cell2mat(par_plan_id);

end