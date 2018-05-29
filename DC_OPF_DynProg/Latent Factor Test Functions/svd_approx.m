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
scen_n = ceil(test_prob.LF_n*test_prob.samp_per)+1;

%% Initialize approximation tools for future use if needed
% mean center scenario data
test_prob.mean_scen_cost = mean(test_prob.scen_op_cost);
A = test_prob.scen_op_cost - test_prob.mean_scen_cost;
%{
% select scenarios from svd latent factors
[~, S, V] = svd(A);
%}
S = test_prob.S;
V = test_prob.V;

best_scen = candexch(V(:,1:LF_n),scen_n, 'tries',10,'display', 'off');

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


residual = std((approx_out(:,best_scen) - cross_val_prob.scen_op_cost(:,best_scen))')';



% replace known values
approx_out(:,best_scen) = cross_val_prob.scen_op_cost(:,best_scen);


tot_out = mean(approx_out,2);
ac_out = mean(cross_val_prob.scen_op_cost,2);
plot((tot_out),'o');
hold on 
plot(ac_out,'x');
for idx = 1:300
    if abs((tot_out(idx) - ac_out(idx))) > (2*residual(idx))
        line([idx idx],[tot_out(idx)+2*residual(idx),tot_out(idx)-2*residual(idx)],'color','r');
    else 
        line([idx idx],[tot_out(idx)+2*residual(idx),tot_out(idx)-2*residual(idx)]); 
    end
end
%plot(mean(cross_val_prob.scen_op_cost,2),'x')



%% Crossvalidate Results
% Compare with pre-calculated crossvalidation set
err_matrix = approx_out - cross_val_prob.scen_op_cost;
err_percent = err_matrix./cross_val_prob.scen_op_cost;

mc_est = mean(cross_val_prob.scen_op_cost(:,randperm(8736,scen_n)),2);

%max err
%best plan difference
%cvar

%problem.second_samp_runtime = cross_val_prob.runtime;

%problem.approx_cost;
%problem.actual_cost;
test_prob.err_out = 100*mean(err_matrix,2)./mean(cross_val_prob.scen_op_cost,2);
test_prob.mean_abs_percent_err = mean(mean(abs(err_percent')));
test_prob.mc_tot_percent_err = mean(abs((mc_est./mean(cross_val_prob.scen_op_cost,2))-1));
test_prob.tot_abs_percent_err = mean(abs((sum(approx_out,2)./(sum(cross_val_prob.scen_op_cost,2)))-1));
test_prob.err_sum_squares = norm(err_matrix, 'fro');
test_prob.r_squared_est = sum(diag(S(1:LF_n,1:LF_n).^2))/sum(diag(S).^2);
test_prob.obs_r_squared = 1 - norm(cross_val_prob.scen_op_cost-approx_out,'fro')/norm(cross_val_prob.scen_op_cost,'fro');
test_prob.obs_r_sq2 = 1 - norm(err_matrix,'fro')^2/norm(cross_val_prob.scen_op_cost-mean(mean(cross_val_prob.scen_op_cost)),'fro').^2;
test_prob.tot_r_sq2 = 1-norm(mean(err_matrix,2))^2/norm(mean(cross_val_prob.scen_op_cost,2)-mean(mean(cross_val_prob.scen_op_cost)))^2;
end