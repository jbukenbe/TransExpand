function [x_fit_id, y_fit] = new_samp_from_beta(problem, beta_op, good_plan_threshold)
% This function takes the problem structure and a beta vector from a
% regression and retuns a vector of plan ids with reasonable estimated fit 

%History            
%Version    Date        Who     Summary
%1          11/18/2017  JesseB  Initial Version
%2          11/30/2017  JesseB  Added Line number limit
%3          11/30/2017  JesseB  Updated rolling best plan matrix
%4          11/30/2017  JesseB  Start with plans with fewer lines

% initialize data
plan_n = problem.params.plan.n;
n_line = problem.params.cand.n;
use_int = problem.params.pls.interaction;
dist_mat_size = problem.params.pls.dist_mat_size;
max_new_lines = problem.params.max_new_lines;
between_cleanup_loop = 12;
x_fit_id = [];
y_fit = [];

% loop through subsections of plans
fit_loop_n = ceil(plan_n / dist_mat_size);
fit_idx = 0;

% search through plans with fewest lines until enough good ones are found
while fit_idx < fit_loop_n
    % initialize fist batch of plans
    x_fit_row = 0;
    k = 0;
    fit_samp_n = problem.params.pls.fit_samp_n;
    while fit_samp_n > x_fit_row
        k = k + 1;
        x_fit_row = x_fit_row + nchoosek(n_line ,k);
        if k == n_line
            fit_samp_n = x_fit_row;
        end
    end
    x_fit = zeros(x_fit_row, n_line)';
    x_fit_full = 0;
    for k_idx = 1:k
        line_idx = nchoosek(1:n_line,k_idx);
        [row_x_fit, col_x_fit] = size(line_idx);
        x_fit_col = n_line*((1:row_x_fit)' + x_fit_full -1);
        x_fit_full = x_fit_full + nchoosek(n_line,k_idx);
        for col_idx = 1:col_x_fit
            x_fit(line_idx(:,col_idx) + x_fit_col)= 1;
        end
    end
    x_fit = [x_fit', zeros(x_fit_row,x_col_n-n_line)];
    
    % initialize parallel data storage cells
    par_x_fit_storage = cell(between_cleanup_loop,1);
    par_y_fit_storage = cell(between_cleanup_loop,1);

    parfor b_idx = 1: between_cleanup_loop
        % find index for full array
        f_idx = fit_idx + b_idx;
        
        % retrive sub set of plans for this iteration
        if f_idx == fit_loop_n
            x_dist_id = (((f_idx-1)*dist_mat_size+1):plan_n)';
        elseif f_idx < fit_loop_n 
            x_dist_id = (((f_idx-1)*dist_mat_size+1):(f_idx*dist_mat_size))';
        else
            x_dist_id = [];
        end
        
        if ~isempty(x_dist_id)
            % convert plan ID to binary lines
            x_fit_dist = de2bi(x_dist_id-1, n_line);

            % remove plans with too many new lines
            x_fit_dist = x_fit_dist(sum(x_fit_dist,2) < max_new_lines,:);

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
            par_x_fit_storage{b_idx} = x_dist_id(logical_good_plan);
            par_y_fit_storage{b_idx} = y_fit_dist(logical_good_plan);
        end  
    end

    % cleanup parallel data
    good_x_id = [x_fit_id; cell2mat(par_x_fit_storage)];
    good_y = [y_fit; cell2mat(par_y_fit_storage)];
    clear par_x_fit_storage par_y_fit_storage

    % remove plans that could not be returned by the function
    [good_y, good_y_id] = sort(good_y);
    x_fit_id = good_x_id(good_y_id);
    samp_size = min(problem.params.pls.fit_samp_n, size(x_fit_id,1));
    x_fit_id = x_fit_id(1: samp_size);
    y_fit = good_y(1: samp_size);
    fit_idx = fit_idx + between_cleanup_loop;
    clear good_y good_y_id good_x_id;
end
end