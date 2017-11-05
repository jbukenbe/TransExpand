function write_struct(problem)
% This function takes the problem structure for the transmission expansion
% problem and saves the relavent output as a matfile. This way results from
% HPC batch runs can be saved and experimented on later

%History            
%Version    Date        Who     Summary
%1          11/02/2017  JesseB  Initial Version

filename = 'output.mat';
output = matfile(filename);
output.plan_id = problem.plan_id;
output.scen_op_cost = problem.scen_op_cost;
output.cand_op_cost = problem.cand_op_cost;
output.cand_full_cost = problem.cand_full_cost;

end