function problem = LF_gen_array(array_id)

%History            
%Version    Date        Who     Summary
%1          07/21/2018  JesseB  Adapted from LF_array_call

%% General Data
% basic data
problem.load_growth = 1.1;
problem.scen_n = 6;

% line data
line_runs = 5;
existing_plan = 1:654;
line_data = matfile('gen_exp_line_samp');
line_list_id = (1:line_runs)';

% gen data
gen_runs = 40;
gen_list_id = ((1:gen_runs)+(gen_runs*(array_id-1)))';
gen_exp_data = matfile('gen_exp_data');
gen_data = matfile('gen_data');
gen_samp = gen_exp_data.gen_samp;
gen_stop = gen_exp_data.gens_built;
gen_cand = gen_exp_data.gen_cand;
gen_built = gen_data.built;
gen_exist = find(gen_built ~=0);


%% This Array Data
array_len = gen_runs*line_runs;
this_array = (1:array_len)+(array_len*(array_id-1));

line_run_list = repmat(line_list_id,gen_runs,1);
gen_run_list = repelem(gen_list_id,line_runs,1);


%% Storage Initialization
outfile_name = sprintf('%s_%d','run',array_id);
output = matfile(outfile_name, 'Writable',true);
output.scen_op_cost = zeros(672,length(this_array));



%% Set Line data for run
problem.scen_w = ones(168,1);
for a_idx = 1:array_len
% find lines that are used in this scenario
    new_lines = double(line_data.plans(line_run_list(a_idx),:)).*(655:705);
    new_lines(new_lines == 0) = [];
    problem.line_plan = [existing_plan, new_lines];   
    
% generators that are used in this scenario
    new_gens = gen_samp(1:gen_stop(gen_run_list(a_idx)),gen_run_list(a_idx));
    new_gens = gen_cand(new_gens);
    problem.gen_plan = sort([gen_exist;new_gens]);
    
    
%% Run Plan          
    for w_idx = 1:4
        scen_offset = (w_idx-1)*(168*12)+1;
        problem.scen_list = scen_offset:(scen_offset+167);
        problem = run_5kgen_gams(problem);   

    % save output for run
        storage_r = (1:168)+(168*(w_idx-1));
        output.scen_op_cost(storage_r,a_idx)= problem.scen_op_cost;
    end
end

end