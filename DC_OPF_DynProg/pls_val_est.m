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
%7          12/01/2017  JesseB  Tweaked to cope with new search parameters
%8          12/02/2017  JesseB  New sample process outsouced to beta functions
%9          12/03/2017  JesseB  Updated to work with new uint8 plan id

%% Initialization
% Extract needed data from problem
z_idx = problem.z_idx;
y = problem.cand_full_cost(problem.samp_range(1,z_idx):problem.samp_range(2,z_idx));
x = problem.samp;
n_line = problem.params.cand.n(z_idx);
use_int = problem.params.pls.interaction;
n_comp = problem.params.pls.n_comp;

% Construct X matrix for regression
if use_int
    x_row_n = size(x,1);
    x_col_n = n_line*(n_line-1)/2;
    x = [x,zeros(x_row_n, x_col_n, 'uint8')];  
end 

% Standardize line data into experimental form -1 = no line 1 = line
%x(:,1:n_line) = (x(:,1:n_line) - .5).*2;

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
[~, ~, ~, ~, beta] = plsregress(double(x),y,n_comp);

%% Output
problem.pls.beta{z_idx} = beta;

end

