function problem = svd_approx(problem, varargin)
% This function takes a generic problem structure input, and outputs either
% the information needed to sample future scenarios and make related
% approximations, or it outputs the approximations to unobserved scenarios
% along with the approximate U values.

%History            
%Version    Date        Who     Summary
%1          02/19/2018  JesseB  Initial Version


%% Initialization Data
lat_fact = problem.lat_fact;

%% Initialize approximation tools for future use if needed
if any(strcmpi('initialize', varargin))
    % mean center scenario data
    problem.mean_scen_cost = mean(problem.scen_op_cost);
    A = problem.scen_op_cost - problem.mean_scen_cost;
    
    % select scenarios from svd latent factors
    [~, S, V] = svd(A);
    [~,factor_id] = sort(abs(V(:,1:lat_fact)),'descend');
    factor_id = unique(factor_id','stable'); 
    best_scen = factor_id(1:lat_fact);

    % find approximation matricies
    X = (S(1:lat_fact,1:lat_fact)*V(best_scen,1:lat_fact)')';
    H = V(:,1:lat_fact)*S(1:lat_fact,1:lat_fact)*inv(X'*X)*X';

    % output data
    problem.params.svd.s_values = S ;
    problem.params.svd.directions = V;
    problem.params.svd.X = X;
    problem.params.svd.hat = H;
    problem.params.svd.latent_scen = best_scen;
    problem.params.svd.filler_scen = setdiff(1:problem.params.scen.n,best_scen);
end

%% Approximate unobserved scenarios if needed
if any(strcmpi('estimate', varargin))
    % mean center partial input vector
    best_scen = problem.params.svd.latent_scen;
    partial_vec = problem.partial_solution;
    partial_vec(:,best_scen) = partial_vec(:,best_scen) - problem.mean_scen_cost(best_scen);
    
    % fill missing values with approximation Hat matrix
    approx_out = problem.params.svd.hat*partial_vec(:,best_scen)' + problem.mean_scen_cost;
    
    % replace known values
    approx_out(:,best_scen) = partial_vec(:,best_scen);
    
    % load output data
    problem.approx_out = approx_out;
end

end