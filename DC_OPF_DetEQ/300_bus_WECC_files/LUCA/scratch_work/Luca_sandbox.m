function problem = Luca_sandbox()
% This demo is a playground function to test the implimentation of LUCA
% as an object oriented program and create a high level program for
% organizing the needed components

%History            
%Version    Date        Who     Summary
%1          11/13/2018  JesseB  Created (not working)
%2          11/19/2018  JesseB  Initial version through 'get gexp sample'
%3          11/27/2018  JesseB  working up to creating tstnep_params 


%% Data
% Load data
load opf_runs_1_to_25 gexp hr_op_cost lexp
hr_n = length(hr_op_cost);
ln_n = max(lexp);
gn_n = max(gexp);
gexp_samp_n = 25;
lf_n = 18;
cluster_n = 11;
lexp_samp_n = 4;
run_n = length(gexp);

luca_data_in = zeros(hr_n, gn_n, ln_n);

for r_idx = 1:run_n
    luca_data_in(:, gexp(r_idx), lexp(r_idx)) = hr_op_cost(:, r_idx); 
end
clear gexp hr_op_cost lexp


%% Initialize Latent Uncertainty Clustering Algorithm (LUCA)
Luca = LucaClass(luca_data_in);


%% Get Optimal (hour, network) Pairs for Generation Expansion Approximation 
[Luca, gexp_samp, Xgn, gsamp_lin_id] = Luca.get_gexp_samp(gexp_samp_n);
load gexp_extrapolate_data.mat gsamp_lin_id gexp_full


Luca.gexp_samp_id = gsamp_lin_id;

Xgn = (Luca.Sgn*(Luca.Vgn'))';
Xgn = Xgn(:,1:Luca.params.dim(2));
Luca.gexp_X = Xgn(gsamp_lin_id,:);

Zgn_mean = mean(Luca.unfold(2, Luca.Z));
Luca.gexp_mean =  Zgn_mean(gsamp_lin_id);

% lookup lines and hours that correspond to gexp_samp
%m = matfile('gexp_extrapolate_data');
%m.gexp_samp = gexp_samp;
%m.hr_samp = gexp_samp(:,1);
%m.tl_samp = gexp_samp(:,2);
%m.Xgn = Xgn;
%m.gsamp_lin_id = gsamp_lin_id;

%load('line_samp')
%load('gen_exp_data')

%m.gexp_full = gexp_full;
%m.line_in = lexp_samp(gexp_samp(:,2), :);


%% Run OPF Optimizations to Approximate Gen Expansion 
% If running in parallel, this is where opf runs would be divided 


%OPF_problem_gexp = OptOPFClass(hr_in, line_in, all_gen);
%output = OPF_problem_gexp.solve();

load gext_runs_1_to_25.mat output;

%% Use OPF Objective to Cluster Gen Expansions
Luca = Luca.get_uhat_gexp(output);

load line_costs.mat line_cost
load tstnep_run_3.mat opt_out
texp_plan = opt_out.line_out;
param_set = cell(5,1);
cluster_n_list = [2 5 9 13 21];
for r_idx = 3:5
    cluster_n = cluster_n_list(r_idx);
    Luca = Luca.cluster_gexp(cluster_n, lf_n);
    
    test = CheckTSTNEPClass(Luca, gexp_full,25);
    test = test.get_opf_runs(texp_plan);
    test = test.check_plan(line_cost);


    %% Get (hour, gexp) Pairs for Each Cluster for the TSTNEP Problem
    %Luca = Luca.get_lexp_samp(lexp_samp_n); 
    %tstnep_params = Luca.make_tstnep_params(gexp_full);                       

    % Save tstnep params for running on cluster
    %param_set{r_idx} = tstnep_params;
end
m = matfile('tstnep_param_data');
m.param_set = param_set;


%% Run TSTNEP Problem
TSTNEP_problem = OptTSTNEPClass(tstnep_params);
TSTNEP_problem = TSTNEP_problem.solve();

% Save results and possibly crossvalidate solution

end

