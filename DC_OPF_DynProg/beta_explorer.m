function [keep_lines, drop_lines] = beta_explorer(beta, line_id)
% This function takes the problem structure and a beta vector from a
% regression returns the lines that are likely to be part of the min cost plans 

%History            
%Version    Date        Who     Summary
%1          12/01/2017  JesseB  Initial Version


%% Initialize data

% make mapping from beta coefficients to lines;
interact_list = nchoosek(line_id,2);
beta_map = [0 0;line_id, zeros(size(line_id,1),1);interact_list];
beta_line_map = [beta, beta_map];
line_beta = beta(2:length(line_id)+1);
%{
for l_idx = 1:length(line_id)
    line_idx = line_id(l_idx);
    top_b_logic = beta_line_map(:,2) == line_idx;
    bot_b_logic = beta_line_map(:,3) == line_idx;
    line_beta(l_idx) = beta_line_map(:,1)'*(top_b_logic+bot_b_logic);
end
%}
% find lines to keep and those to drop
keep_lines = line_find(1, line_beta, line_id, beta_line_map);
drop_lines = line_find(0, line_beta, line_id, beta_line_map);

%% Run the Algorithm
function return_lines = line_find(keep, line_beta, line_id, beta_line_map)
    % find line with most impactful beta coefficient depending on if we are
    % looking for lines to keep or drop
    if keep
        [this_beta ,return_line] = min(line_beta);
    else
        [this_beta ,return_line] = max(line_beta);
    end

    % store this line in the list to be returned
    return_line_id = line_id(return_line);
    return_lines = return_line_id;
    stop = 0;
    
    % loop through algorithm until no new lines are chosen
    while ~stop
        % add interaction beta to line beta based on lines already selected
        beta_up_top = beta_line_map(beta_line_map(:,3) == return_line_id,1);
        beta_up_bot = beta_line_map(beta_line_map(:,2) == return_line_id,1);
        
        % prevent interaction factors from being double counted
        beta_line_map(beta_line_map(:,2) == return_line_id,1) = 0;
        beta_line_map(beta_line_map(:,3) == return_line_id,1) = 0;

        % update beta function and find next line
        if keep 
            beta_up = [beta_up_top; -this_beta; beta_up_bot(2:end)];
            line_beta = line_beta + beta_up;
            [this_beta ,return_line] = min(line_beta);
        else
            beta_up = [beta_up_top; this_beta; beta_up_bot(2:end)];
            line_beta = line_beta - beta_up;
            [this_beta ,return_line] = max(line_beta);
        end
        
        % see if next line is new and stop loop if not
        return_line_id = line_id(return_line);
        stop = sum(return_lines == return_line_id);

        if ~stop
            return_lines = [return_lines;return_line_id];
        end
    end
end
end
