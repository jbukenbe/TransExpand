function [unique_opt_plans,map_to_unique , map_to_original_plan] = array_extract()
% This function extracts the optimal plan line data from the hpc test array
% and saves a matfile with the index of candidate lines for each unique
% plan. The map back to the original plan order and into the unique set are
% also returned

%History            
%Version    Date        Who     Summary
%1          06/11/2018  JesseB  Initial Version

opt_plans = zeros(320,51);
fill_idx = 0;
for a_idx = 2:6
    outfile_name = sprintf('%s_%d','array',a_idx);
    load(outfile_name);
        
    for p_idx = 1:length(lines_built)
        opt_plans(p_idx+fill_idx,lines_built{p_idx}) = 1;
    end
    fill_idx = fill_idx + p_idx;
end
[unique_opt_plans,map_to_unique,map_to_original_plan] = unique(opt_plans,'rows');
n = size(unique_opt_plans,1);
out_name = sprintf('%ss_1_to_%d','opt_plan',n);
m = matfile(out_name);
m.plans = unique_opt_plans;
m.map_to_original_plan = map_to_original_plan;

end
