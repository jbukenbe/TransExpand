function [x_fit_id, y_fit] = new_samp_from_beta(problem)
% This function takes the problem structure and a beta vector from a
% regression and retuns a vector of plan ids with reasonable estimated fit 

%History            
%Version    Date        Who     Summary
%1          11/18/2017  JesseB  Initial Version
%2          11/30/2017  JesseB  Added Line number limit
%3          11/30/2017  JesseB  Updated rolling best plan matrix
%4          11/30/2017  JesseB  Start with plans with fewer lines
%5          12/01/2017  JesseB  Revamped to pick best lines first
%6          12/02/2017  JesseB  Returns plans with lowest line costs

% initialize data
beta = problem.pls.beta;
good_plan_threshold = problem.pls.good_plan_threshold;
line_id = problem.params.cand.line_id;
n_line = problem.params.cand.n;
use_int = problem.params.pls.interaction;
dist_mat_size = problem.params.pls.dist_mat_size;
max_new_lines = problem.params.max_new_lines;
line_cost = problem.params.new_line_cost;
line_budget = problem.params.line_budget;
between_cleanup_loop = 12;
x_fit_id = [];
line_cost_list = [];
y_fit = [];

% sort best lines by beta vector
best_lines = beta_graph_search(problem);
%best_lines = beta_explorer(problem);

% subset of best lines to look through
sub_best_lines_n = min(length(best_lines),max_new_lines);
subset_best_lines_id = best_lines(1:sub_best_lines_n); 
sub_plan_n = 2^sub_best_lines_n;

% subset best lines index in original problem structure
[~, subset_best_lines] = ismember(subset_best_lines_id, line_id);

% loop through subsections of plans
fit_loop_n = ceil(sub_plan_n / dist_mat_size);
fit_idx = 0;

% search through plans with fewest lines until enough good ones are found
while fit_idx < fit_loop_n    
    % initialize parallel data storage cells
    par_x_fit_storage = cell(between_cleanup_loop,1);
    par_line_cost_storage =cell(between_cleanup_loop,1);
    par_y_fit_storage = cell(between_cleanup_loop,1);

    parfor b_idx = 1: between_cleanup_loop
        % find index for full array
        f_idx = fit_idx + b_idx;
        
        % retrive sub set of plans for this iteration
        if f_idx == fit_loop_n
            x_subset_id = (((f_idx-1)*dist_mat_size+1):sub_plan_n)';
        elseif f_idx < fit_loop_n 
            x_subset_id = (((f_idx-1)*dist_mat_size+1):(f_idx*dist_mat_size))';
        else
            x_subset_id = [];
        end
        
        if ~isempty(x_subset_id)
            % convert plan ID to binary lines
            x_subset = de2bi(x_subset_id-1, sub_best_lines_n);

            % remove plans with too many new lines 
            x_subset = x_subset(sum(x_subset,2) < max_new_lines,:);          
            
            % convert plans from best line form to original problem binary
            x_fit = zeros(size(x_subset,1), n_line);
            x_fit(:,subset_best_lines) = x_subset; 
            
            % remove plans over budget
            line_cost_distributed = x_fit*line_cost;
            x_fit(line_cost_distributed > line_budget,:) = [];
            line_cost_distributed(line_cost_distributed > line_budget,:) = [];
            
            %store plan id in original problem format
            x_original_id = bi2de(x_fit)+1;
            
            % convert plans to experimental form
            %x_fit(:,1:sub_best_lines_n) = (x_fit(:,1:sub_best_lines_n) - .5)*2;
            if use_int
                x_col_start = n_line + 1;    
                for l_idx = 1:(n_line-1)    
                    x_col_end = x_col_start + ((n_line-1) - l_idx);        
                    x_fit(:, x_col_start:x_col_end) = repmat(x_fit(:,l_idx),1,(1+x_col_end-x_col_start)).*x_fit(:,(l_idx+1):n_line);        
                    x_col_start = x_col_end + 1;        
                end    
            end

            % check fit for plans
            y_fit_distributed = [ones(size(x_fit,1),1), x_fit]*beta;
            logical_good_plan = (y_fit_distributed < good_plan_threshold);

            % store plans with best fits
            par_x_fit_storage{b_idx} = x_original_id(logical_good_plan);
            par_line_cost_storage{b_idx} = line_cost_distributed(logical_good_plan);
            par_y_fit_storage{b_idx} = y_fit_distributed(logical_good_plan);
        end  
    end

    % cleanup parallel data
    good_x_id = [x_fit_id; cell2mat(par_x_fit_storage)];
    good_line_cost = [line_cost_list; cell2mat(par_line_cost_storage)];
    good_y = [y_fit; cell2mat(par_y_fit_storage)];
    clear par_x_fit_storage par_y_fit_storage

    % remove plans that could not be returned by the function
    %[good_y, good_y_id] = sort(good_y);
    [good_line_cost, good_line_cost_id] = sort(good_line_cost);
    x_fit_id = good_x_id(good_line_cost_id);
    good_y = good_y(good_line_cost_id);
    samp_size = min(problem.params.pls.fit_samp_n, size(x_fit_id,1));
    x_fit_id = x_fit_id(1: samp_size);
    line_cost_list = good_line_cost(1:samp_size);
    y_fit = good_y(1: samp_size);
    fit_idx = fit_idx + between_cleanup_loop;
    clear good_y good_y_id good_x_id;
end
end
