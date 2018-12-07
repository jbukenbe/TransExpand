function [unique_opt_plans,map_to_unique,map_to_original_plan] = array_extract(exp_string ,line_cost)
% This function extracts the optimal plan line data from the hpc test array
% and saves a matfile with the index of candidate lines for each unique
% plan. The map back to the original plan order and into the unique set are
% also returned

%History            
%Version    Date        Who     Summary
%1          06/11/2018  JesseB  Initial Version


run_n = 32*8;

opt_plans = zeros(run_n,51);
obj_log = zeros(run_n,1);
opt_log = zeros(run_n,1);
time_log = zeros(run_n,1);
bound_log = zeros(run_n,1);

run_list = [1, 3:6, 8:14, 16, 18:22, 24:25, 27, 29:37, 39:40];

fill_idx = 0;
for r_idx = 1:32
    a_idx = run_list(r_idx);
    outfile_name = sprintf('%s_%d',exp_string,a_idx);
    load(outfile_name);
    
    runs_in_array = length(lines_built);
    log_idx =(1:runs_in_array)+(runs_in_array*(a_idx-1));
    
    obj_log(log_idx) = obj_val;
    opt_log(log_idx) = opt_gap;
    time_log(log_idx) = run_time;
    bound_log(log_idx) = lower_bound;
    
    
    for p_idx = 1:runs_in_array
        opt_plans(p_idx+fill_idx,lines_built{p_idx}) = 1;
    end
    fill_idx = fill_idx + p_idx;
end
[unique_opt_plans,map_to_unique,map_to_original_plan] = unique(opt_plans,'rows');
n = size(unique_opt_plans,1);
out_name = sprintf('%s_1_to_%d',exp_string,run_n);
m = matfile(out_name);

m.plans = unique_opt_plans;
m.map_to_original_plan = map_to_original_plan;

m.scen_n = repelem([2;4;8;12;24],20,1);
m.plan_cost = opt_plans*line_cost;
m.opt_plans = opt_plans;
m.obj_val = obj_log;
m.opt_gap = opt_log;
m.run_time = time_log;
m.lower_bound = bound_log;

end
