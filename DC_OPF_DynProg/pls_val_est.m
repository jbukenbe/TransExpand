function problem = pls_val_est(problem)
% This function takes the problem from the dyn_DC_OPF function and uses
% Partial Least Squares (also called Projected on Latent Structure) Regression
% to make estimates of the value of future plans 

%History            
%Version    Date        Who     Summary
%1          10/06/2017  JesseB  Initial Version

%% Initialization
% Extract needed data from problem
ordering = randperm(1024)';
y = problem.cand_cost(ordering(1:50));
%c_line = problem.line_cost;
x_id = problem.plan_id(ordering(1:50));
n_line = problem.params.cand.n;
use_int = problem.interaction;
n_comp = 10;
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
[xl, yl, xs, ys, beta, pctvar] = plsregress(x,y,n_comp);
plot(1:10,cumsum(100*pctvar(2,:)),'-bo')
%% Post Processing
% Extract n_pls_plan new plans from regression and n_rand_plan plans randomly 



%% Output




end