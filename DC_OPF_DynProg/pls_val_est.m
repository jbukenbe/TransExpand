function problem = pls_val_est(problem)
% This function takes the problem from the dyn_DC_OPF function and uses
% Partial Least Squares (also called Projected on Latent Structure) Regression
% to make estimates of the value of future plans 

%History            
%Version    Date        Who     Summary
%1          10/06/2017  JesseB  Initial Version
%2          10/07/2017  JesseB  Added simple search for new plans
%3          10/08/2017  JesseB  Added line costs
%4          11/04/2017  JesseB  Reworked searching for finding min cost plan
%5          11/15/2017  JesseB  Modified to output refined sample for real time search
%6          11/18/2017  JesseB  Pulled out good fit id picking process for large arrays

% To Do:    Avoid duplicate plans

%% Initialization
% Extract needed data from problem
y = problem.cand_op_cost;
line_cost = problem.params.new_line_cost;
x_id = problem.plan_id;
n_line = problem.params.cand.n;
n_plans = 2^n_line;
use_int = problem.params.pls.interaction;
n_comp = problem.params.pls.n_comp;

% Construct X matrix for regression
x_row_n = size(x_id,1);
x_col_n = n_line;
if use_int
    x_col_n = x_col_n*(x_col_n+1)/2;
end 
x = zeros(x_row_n, x_col_n);
x(:,1:n_line) = de2bi(x_id-1, n_line);

% Standardize line data into experimental form -1 = no line 1 = line
x(:,1:n_line) = (x(:,1:n_line) - .5).*2;

% Include first order interactions in X matrix
if use_int
    x_col_start = n_line + 1;
    for l_idx = 1:(n_line-1)        
        x_col_end = x_col_start + ((n_line-1) - l_idx);
        x(:,x_col_start:x_col_end) = repmat(x(:,l_idx),1,(1+x_col_end-x_col_start)).*x(:,(l_idx+1):n_line);
        x_col_start = x_col_end + 1;
    end
end

%% Run PLS regression
[~, ~, ~, ~, beta_op, pctvar] = plsregress(x,y,n_comp);


%% Post Processing
good_plan_threshold = min(y)*1.05;
[x_fit_id, y_fit] = new_samp_from_beta(problem, beta_op, good_plan_threshold);
x_fit_line = de2bi(x_fit_id-1, n_line);

% Get fits for all plans operating costs
[~,fit_rank]=sort(y_fit);
plot(y_fit(fit_rank))

% Isolate plans with lowest estimated operating costs
best_op_set = fit_rank(1:problem.params.pls.line_samp_n);

% Calculate line costs for best plans
plan_line_cost = x_fit_line(best_op_set,:)*line_cost;

% Sort by line cost
[~, best_line_subset] = sort(plan_line_cost);
best_op_set = best_op_set(best_line_subset);

% Calculate full cost for best sets
best_full_set = best_op_set(1:problem.params.pls.refine_samp_n);


%% Output
problem.pls.beta = beta_op;
problem.samp_id = bi2de(x_fit_line(best_full_set,:))+1;


%% Plotting in debug
%{
plot(1:n_comp,cumsum(100*pctvar(2,:)),'-bo')

x_fit_id = problem.plan_id;
x_fit = de2bi(x_fit_id-1, n_line);
x_fit = (x_fit - .5).*2;
if use_int
    for l_idx = 1:(n_line-1)
        x_fit = horzcat(x_fit,x_fit(:,l_idx).*x_fit(:,(l_idx+1):n_line));    
    end
end

y_fit = [ones(size(x_fit,1),1), x_fit]*beta_op+de2bi(x_fit_id-1, n_line)*c_line;
figure
scatter(y_fit, problem.cand_op_cost+de2bi(x_fit_id-1, n_line)*c_line);
xlabel('fitted value'); ylabel('Actual Plan Cost'); title('Example Predictive Power from 300 Sample Plans'); 
%}




end

