function problem = LF_array_call(array_id)
% This function is called by the HPC array. The array id is used to
% construct the problem parameters for the current experimental settup

%History            
%Version    Date        Who     Summary
%1          04/03/2018  JesseB  Initial Version
%2          05/22/2018  JesseB  Broke for testing new GAMS model feeders
%3          07/19/2018  JesseB  Adapted to LF test array solutions

%% User Defined Settings


%% Model Settings
% basic settings
problem.week = 1;
problem.bus_n = 312;
problem.gen_n = 980;
problem.non_dispatch_gen = 332:980;
problem.load_growth = 1.3;

% candidate plan line settings
problem.existing_plan = 1:654;
line_data = matfile('other_234.mat');
is_line_data = matfile('is_line_data');
run_list = reshape(1:234,13,18);
run_list = run_list(:,14);
[line_runs,~] = size(run_list);
line_runs = 1;



% scenario settings
problem.scen_n = 6;

%% Storage Initialization
outfile_name = sprintf('%s_%d','last_run',array_id);
output = matfile(outfile_name, 'Writable',true);
output.scen_op_cost = zeros(8736,line_runs);
output.shed_gen = zeros(8736,line_runs);
output.pns = zeros(8736,line_runs);
output.lines = zeros(line_runs,51);


%% Set Line data for run
for l_idx = 1:line_runs
%    problem = run_bender_tep_gams(problem);
    %problem.new_lines = problem.lines+654;
    if array_id <5
        problem.new_lines = double(is_line_data.plans(array_id,:)).*(655:705);
        problem.new_lines(problem.new_lines == 0) = [];
        problem.candidate_plan = [problem.existing_plan,problem.new_lines];        
        output.lines(l_idx,1:51) = is_line_data.plans(array_id,:);
    
    else
        run_idx = run_list(array_id-4);
        problem.new_lines = double(line_data.plans(run_idx,:)).*(655:705);
        problem.new_lines(problem.new_lines == 0) = [];
        problem.candidate_plan = [problem.existing_plan,problem.new_lines];        
        output.lines(l_idx,1:51) = line_data.plans(run_idx,:);
    end
%% Run Plan          
    for w_idx = 1:52
        problem.scen_offset = (w_idx-1)*168+1;
        problem = run_opf_gams(problem);   

    % save output for run
        storage_r = (168*(w_idx-1)+1):168*w_idx;
        output.scen_op_cost(storage_r,l_idx)= problem.scen_op_cost;
        output.shed_gen(storage_r,l_idx)= problem.shed_gen;
        output.pns(storage_r,l_idx)= problem.pns;
        
    end
end

end