function problem = par_plan_ops(problem)
% This function takes the problem structure for the transmission expansion
% problem and solves the plan costs in parallel 

%History            
%Version    Date        Who     Summary
%1          10/07/2017  JesseB  Adapted from dyn_pls_demo
%2          12/03/2017  JesseB  Minor restructuring to work with full sample as unit8
%3          12/03/2017  JesseB  Added SVD approximator first version
%4          02/19/2017  JesseB  Compatible with integer lines

%% Load needed problem parameters
params = problem.params;
samp = problem.samp;
z_idx = problem.z_idx;
samp_size = size(samp,1);
line_id = params.cand.line_id{z_idx};
dec_built = params.line.dec_built{z_idx};

%% Initialize Storage
par_scen_op_cost = cell(samp_size, params.scen.n);
par_cand_op_cost = cell(samp_size, 1);
par_cand_full_cost = cell(samp_size,1);
    
%% Run all sample plans in parallel
parfor c_idx = 1:samp_size
% load relavent data for DC OPF for this plan
    new_line_idx = line_id.*double(samp(c_idx,:)');
    new_line_idx = nonzeros(new_line_idx);
    dec_lines = dec_built;
    dec_lines(new_line_idx) = true;
    new_number_built = zeros(size(dec_lines,1),1);
    new_number_built(new_line_idx) = nonzeros(double(samp(c_idx,:)'));
    cand_line_n = params.line.per_corridor + new_number_built;

% copy parameter data for parallel running
    par_params = params;
    par_params.cand.imp = params.line.imp(dec_lines);
    par_params.cand.res = params.line.res(dec_lines);
    par_params.cand.from_to = params.line.from_to(dec_lines,:);
    par_params.cand.per_corridor = cand_line_n(dec_lines,:);
    par_params.cand.max_flow = params.line.max_flow(dec_lines);

    total_line_cost = sum(par_params.line.cost(new_line_idx));
    
%% Run Scenario Optimization
% run optimization to get operating cost of this plan for each scenario
    op_cost = cell(1,par_params.scen.n);
    op_cost(1,:) = {0};
    if ~par_params.svd.use_latent_fac
        % calculate all operating scenarios
        for scen = 1:par_params.scen.n
            op_cost{scen} = DC_OPF_operations(par_params, dec_lines, scen);
        end
    else
        % estimate operations with svd latent factor approximation
         for scen = 1:par_params.svd.scen_n
            scen_latent = par_params.svd.latent_scen(scen);
            op_cost{scen_latent} = DC_OPF_operations(par_params, dec_lines, scen_latent);
         end
         
        % approximate total cost
        approx_out = mat_fill_svd(par_params.svd.s_values, par_params.svd.directions, cell2mat(op_cost));
        approx_cost = approx_out*par_params.scen.p;
        
        % if this looks like a good plan, calculate the full cost
        if approx_cost < par_params.svd.good_plan_threshold
            for scen = 1:(par_params.scen.n - par_params.svd.scen_n)
                scen_filler = par_params.svd.filler_scen(scen);
                op_cost{scen_filler} = DC_OPF_operations(par_params, dec_lines, scen_filler);
            end
        else
            op_cost = num2cell(approx_out);
        end
    end
    
    
%% Store Output Data
    par_scen_op_cost(c_idx,:) = op_cost;
    op_cost = cell2mat(op_cost);
    par_cand_op_cost{c_idx} = op_cost*par_params.scen.p;
    par_cand_full_cost{c_idx} = par_cand_op_cost{c_idx} + total_line_cost;
end
        
%% Parallel Data Cleanup
problem.scen_op_cost = cell2mat(par_scen_op_cost);
problem.cand_op_cost = cell2mat(par_cand_op_cost);
problem.cand_full_cost = cell2mat(par_cand_full_cost); 

end