function problem = LF_test_init(init_samp)
% This fuction runst the initial trials for the SVD latent factor
% approximator algorithm

%History            
%Version    Date        Who     Summary
%1          02/19/2018  JesseB  Initial Version

%% Initialization
% load problem data for LF settings
problem.z_idx = 1; problem.problem_size = 196;
problem.params = DC_OPF_init(problem);    
problem.params.initial_samp_n = init_samp;

%% Initial Plan Sample
%problem.samp_id = randperm(plan_n,params.initial_samp_n)';
problem.samp = fraction_fact_samp(problem);
problem.plan{1} = problem.samp;
problem.samp_range(1,problem.z_idx) = 1;
problem.samp_range(2,problem.z_idx) = size(problem.samp,1); 

%% Run OPF
run_time = tic;
problem = par_plan_ops(problem);
problem.run_time = toc(run_time);

end