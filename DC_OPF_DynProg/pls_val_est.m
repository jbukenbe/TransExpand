function problem = pls_val_est(problem)
% This function takes the problem from the dyn_DC_OPF function and uses
% Partial Least Squares (also called Projected on Latent Structure) Regression
% to make estimates of the value of future plans 

%History            
%Version    Date        Who     Summary
%1          10/06/2017  JesseB  Initial Version
%2          10/07/2017  JesseB  Added simple search for new plans

% To Do:    Include Line Costs
%           Avoid duplicate plans

%% Initialization
% Extract needed data from problem
y = problem.cand_cost;
%c_line = problem.line_cost;
x_id = problem.plan_id;
n_line = problem.params.cand.n;
use_int = problem.params.interaction;
n_comp = problem.params.n_comp;
% Construct X matrix for regression
x = de2bi(x_id-1, n_line);

% Standardize line data into experimental form -1 = no line 1 = line
x = (x - .5).*2;

% Include first order interactions in X matrix
if use_int
    for l_idx = 1:(n_line-1)
        x = horzcat(x,x(:,l_idx).*x(:,(l_idx+1):n_line));    
    end
end
%% Run PLS regression
[~, ~, ~, ~, beta, pctvar] = plsregress(x,y,n_comp);


%% Post Processing
% Extract inform_sample_n new plans from regression and explore_sample_n plans randomly 
n_pls_plan = problem.params.inform_sample_n;
good_lines = beta(2:n_line+1)' < 0;
good_lines = (good_lines-.5).*2;
new_inform_plans = good_lines(ones(n_pls_plan,1),:);
[~,line_rank]=sort(beta(2:(n_line+1)).^2);
plan_flip = -de2bi(0:(n_pls_plan-1),n_line);
plan_flip = plan_flip(:,line_rank);
plan_flip = (plan_flip+.5).*2;

new_inform_plans = bi2de((new_inform_plans.*plan_flip+1)/2)+1;

new_explore_plans = randperm(2^n_line,problem.params.explore_sample_n)';
%% Output
problem.new_plan_id = [new_inform_plans; new_explore_plans];



%% Plotting in debug

plot(1:n_comp,cumsum(100*pctvar(2,:)),'-bo')
%{
x_fit_id = problem.plan_id;
x_fit = de2bi(x_fit_id-1, n_line);
x_fit = (x_fit - .5).*2;
if use_int
    for l_idx = 1:(n_line-1)
        x_fit = horzcat(x_fit,x_fit(:,l_idx).*x_fit(:,(l_idx+1):n_line));    
    end
end

y_fit = [ones(size(x_fit,1),1), x_fit]*beta;
figure
scatter(y_fit, problem.cand_cost);
xlabel('fitted value'); ylabel('Actual Plan Cost'); title('Example Predictive Power from 300 Sample Plans'); 
%}
end

