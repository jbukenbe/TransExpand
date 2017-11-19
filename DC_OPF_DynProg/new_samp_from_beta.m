function [x_fit_id, y_fit] = new_samp_from_beta(problem, beta_op, good_plan_threshold)
% This function takes the problem structure and a beta vector from a
% regression and retuns a vector of plan ids with reasonable estimated fit 

%History            
%Version    Date        Who     Summary
%1          10/06/2017  JesseB  Initial Version


% initialize data
plan_n = problem.params.plan.n;
n_line = problem.params.cand.n;
use_int = problem.params.pls.interaction;
dist_mat_size = problem.params.pls.dist_mat_size;

% loop through subsections of plans
fit_loop_n = ceil(plan_n / dist_mat_size);
par_x_fit_storage = cell(fit_loop_n,1);
par_y_fit_storage = cell(fit_loop_n,1);

parfor f_idx = 1: fit_loop_n
    
% retrive sub set of plans for this iteration
    if f_idx == fit_loop_n
        x_dist_id = (((f_idx-1)*dist_mat_size+1):plan_n)';
    else
        x_dist_id = (((f_idx-1)*dist_mat_size+1):(f_idx*dist_mat_size))';
    end
x_fit_dist = de2bi(x_dist_id-1, n_line);
    
% convert subset to experimental form
x_fit_dist(:,1:n_line) = (x_fit_dist(:,1:n_line) - .5)*2;
if use_int
    x_col_start = n_line + 1;
    for l_idx = 1:(n_line-1)
        x_col_end = x_col_start + ((n_line-1) - l_idx);
        x_fit_dist(:, x_col_start:x_col_end) = repmat(x_fit_dist(:,l_idx),1,(1+x_col_end-x_col_start)).*x_fit_dist(:,(l_idx+1):n_line);
        x_col_start = x_col_end + 1;
    end
end

% check fit for subset of plans
y_fit_dist = [ones(size(x_fit_dist,1),1), x_fit_dist]*beta_op;
logical_good_plan = (y_fit_dist < good_plan_threshold);

% store plans with best fits
par_x_fit_storage{f_idx} = x_dist_id(logical_good_plan);
par_y_fit_storage{f_idx} = y_fit_dist(logical_good_plan);
end

% cleanup parallel data
good_x_id = cell2mat(par_x_fit_storage);
good_y = cell2mat(par_y_fit_storage);

% output best plan ids
[good_y, good_y_id] = sort(good_y);
x_fit_id = good_x_id(good_y_id);
samp_size = min(problem.params.pls.fit_samp_n, size(x_fit_id,1));
x_fit_id = x_fit_id(1: samp_size);
y_fit = good_y(1: samp_size);

end