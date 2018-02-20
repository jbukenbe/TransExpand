function test_prob = svd_approx(test_prob, cross_val_prob)
% This function takes a generic problem structure input, and outputs either
% the information needed to sample future scenarios and make related
% approximations, or it outputs the approximations to unobserved scenarios
% along with the approximate U values.

%History            
%Version    Date        Who     Summary
%1          02/19/2018  JesseB  Initial Version


%% Initialization Data
LF_n = test_prob.LF_n;
scen_n = test_prob.LF_n*test_prob.samp_per;

%% Initialize approximation tools for future use if needed
% mean center scenario data
test_prob.mean_scen_cost = mean(test_prob.scen_op_cost);
A = test_prob.scen_op_cost - test_prob.mean_scen_cost;

% select scenarios from svd latent factors
[~, S, V] = svd(A);
[~,factor_id] = sort(abs(V(:,1:LF_n)),'descend');
factor_id = unique(factor_id','stable'); 
best_scen = factor_id(1:scen_n);

% find approximation matricies
X = (S(1:LF_n,1:LF_n)*V(best_scen,1:LF_n)')';
%H = V(:,1:LF_n)*S(1:LF_n,1:LF_n)*inv(X'*X)*X';
H = V(:,1:LF_n)*S(1:LF_n,1:LF_n)*((X'*X)\(X'));

% output data
test_prob.svd.s_values = S ;
test_prob.svd.directions = V;
test_prob.svd.X = X;
test_prob.svd.hat = H;
test_prob.svd.latent_scen = best_scen;
test_prob.svd.filler_scen = setdiff(1:test_prob.params.scen.n,test_prob.svd.latent_scen);


%% Approximate unobserved scenarios if needed
% create simulated secondary set from crossval problem
partial_vec = cross_val_prob.scen_op_cost;
partial_vec(:,test_prob.svd.filler_scen) = 0;

% mean center partial input vector
partial_vec(:,best_scen) = partial_vec(:,best_scen) - test_prob.mean_scen_cost(best_scen);

% fill missing values with approximation Hat matrix
approx_out = (H*partial_vec(:,best_scen)')' + test_prob.mean_scen_cost;

% replace known values
approx_out(:,best_scen) = cross_val_prob.scen_op_cost(:,best_scen);


%% Crossvalidate Results
% Compare with pre-calculated crossvalidation set
err_matrix = approx_out - cross_val_prob.scen_op_cost;
err_percent = err_matrix./cross_val_prob.scen_op_cost;

%max err
%best plan difference
%cvar

%problem.second_samp_runtime = cross_val_prob.runtime;

%problem.approx_cost;
%problem.actual_cost;
test_prob.mean_abs_percent_err = mean(max(abs(err_percent')));
test_prob.sum_squares_est = norm(err_matrix, 'fro');
test_prob.r_squared_est = sum(diag(S(1:LF_n,1:LF_n)))/sum(diag(S));
test_prob.obs_r_squared =1-(test_prob.err_sum_squares/norm(cross_val_prob.scen_op_cost,'fro'));

end