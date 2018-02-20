function problem = LF_crossval_set(problem,varargin)
% This function will create the crossvalidation data for the LF
% approximation algorithm

%History            
%Version    Date        Who     Summary
%1          02/19/2018  JesseB  Initial Version

%% Calculate Data if needed
if any(strcmpi('calculate', varargin))

problem.z_idx = 1; problem.problem_size = 15;
problem.params = DC_OPF_init(problem);    
problem.params.initial_samp_n = problem.samp_n;

problem.samp = fraction_fact_samp(problem);
problem.plan{1} = problem.samp;
problem.samp_range(1,problem.z_idx) = 1;
problem.samp_range(2,problem.z_idx) = size(problem.samp,1); 

% run opf
problem = par_plan_ops(problem);
end


%% Load Data if available
if any(strcmpi('load', varargin))


end

end