function [best_plan_full_val, best_plan_id] = pls_val_est(problem)
% This function takes the problem from the dyn_DC_OPF function and uses
% Partial Least Squares (also called Projected on Latent Structure) Regression
% to make estimates of the value of future plans 

%History            
%Version    Date        Who     Summary
%1          10/06/2017  JesseB  Initial Version
%2          10/07/2017  JesseB  Added simple search for new plans
%3          10/08/2017  JesseB  Added line costs
%4          11/04/2017  JesseB  Reworked searching for finding min cost plan

% To Do:    Include Line Costs
%           Avoid duplicate plans

%% Initialization
% Extract needed data from problem
y = problem.cand_op_cost;
line_cost = problem.params.new_line_cost;
x_id = problem.plan_id;
n_line = problem.params.cand.n;
n_plans = 2^n_line;
use_int = problem.params.interaction;
n_comp = problem.params.n_comp;

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
        x(:,x_col_start:x_col_end) = x(:,l_idx).*x(:,(l_idx+1):n_line);
        x_col_start = x_col_end + 1;
    end
end

%% Run PLS regression
[~, ~, ~, ~, beta_op, pctvar] = plsregress(x,y,n_comp);


%% Post Processing
% create full matrix of plans for fit estimation
x_fit = zeros(n_plans, x_col_n);
x_fit(:,1:n_line) = de2bi(0:(n_plans-1), n_line);
x_fit(:,1:n_line) = (x_fit(:,1:n_line) - .5).*2;
if use_int
    x_col_start = n_line + 1;
    for l_idx = 1:(n_line-1)
        x_col_end = x_col_start + ((n_line-1) - l_idx);
        x_fit(:, x_col_start:x_col_end) = x_fit(:,l_idx).*x_fit(:,(l_idx+1):n_line);
        x_col_start = x_col_end + 1;
    end
end

% Get fits for all plans operating costs
y_fit = [ones(size(x_fit,1),1), x_fit]*beta_op;
[~,fit_rank]=sort(y_fit);
%plot(problem.real_op(fit_rank(1:1000000)))

% Isolate best plans
best_op_set = fit_rank(1:100000);

% Calculate line costs for best plans
plan_line_cost = x_fit(best_op_set,1:n_line)*line_cost;

% Sort by line cost
[~, best_line_subset] = sort(plan_line_cost);
best_op_set = best_op_set(best_line_subset);

% Calculate full cost for best sets
best_full_set = best_op_set(1:500);
[best_full_val, best_full_subset] = sort(problem.cand_full_cost(best_full_set));
%scatter(1:size(best_full_val,1),best_full_val)

%% Output
best_plan_id = best_full_set(best_full_subset(1:100));
best_plan_full_val = best_full_val(1:100); 

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

