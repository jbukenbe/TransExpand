function problem = LF_array_call(array_id)
% This function is called by the HPC array. The array id is used to
% construct the problem parameters for the current experimental settup

%History            
%Version    Date        Who     Summary
%1          04/03/2018  JesseB  Initial Version

%% User Defined Settings

tic
%% Model Settings
% basic settings
problem.week = 1;
problem.bus_n = 312;
problem.gen_n = 980;
problem.non_dispatch_gen = 332:980;
problem.load_growth = 1.3;

% candidate plan line settings
line_runs = 12;
problem.existing_plan = 1:654;
line_data = matfile('line_samp_two.mat');
run_list = reshape(1:300,12,25);

% scenario settings
problem.scen_n = 6;

%% Storage Initialization
outfile_name = sprintf('%s_%d','line',array_id);
output = matfile(outfile_name, 'Writable',true);
output.scen_op_cost = zeros(8736,line_runs);

%% Set Line data for run
for l_idx = 1:line_runs
    problem.new_lines = double(line_data.samp(run_list(l_idx,array_id),:)).*(655:705);
    problem.new_lines(problem.new_lines == 0) = [];
    problem.candidate_plan = [problem.existing_plan,problem.new_lines];        
    
%% Run Plan
    for w_idx = 1:52
        problem.scen_offset = (w_idx-1)*168+1;
        problem = run_plan_in_gams(problem);

    % save output for run
        storage_r = (168*(w_idx-1)+1):168*w_idx;
        output.scen_op_cost(storage_r,l_idx)= problem.scen_op_cost;
    end
end
runtime = toc
end