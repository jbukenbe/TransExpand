function problem = new_samp_from_beta(problem)
% This function takes the problem structure and a beta vector from a
% regression and retuns a vector of plan ids with reasonable estimated fit 

%History            
%Version    Date        Who     Summary
%1          11/18/2017  JesseB  Initial Version

%% Initialize Data
z_idx = problem.z_idx;
line_id = problem.params.cand.line_id{z_idx-1};
beta = problem.pls.beta{z_idx-1};
full_line_id = problem.params.cand.line_id{1};
lock_on = problem.lock_on{z_idx-1};
lock_off = problem.lock_off{z_idx-1};

lock_off_n = length(lock_off);
lock_off(randperm(lock_off_n, ceil(lock_off_n*.10))) = [];

%% Find Lines with Good Beta
if problem.params.pls.interaction 
    %best_lines = beta_graph_search(problem);
    [keep_lines, drop_lines] = beta_explorer(beta,line_id);
else
    [~,beta_order] = sort(beta(2:end));
    keep_lines = line_id(beta_order);
    drop_lines = line_id(flip(beta_order));
end
    
%% Select Lines to Lock and Sample
drop_n = length(drop_lines);
lock_off = [lock_off;drop_lines(1:(ceil(.5*drop_n)))];
new_line_id = setdiff(full_line_id, lock_off);
        
%% Update Problem Structure
problem.params.cand.line_id{z_idx} = new_line_id;
problem.params.new_line_cost{z_idx} = problem.params.line.cost(new_line_id);
problem.params.cand.n(z_idx) = length(new_line_id);
logic_lock_on = zeros(problem.params.line.n,1,'logical');
logic_lock_on(lock_on) = 1;
problem.params.line.dec_built{z_idx} = logical(logic_lock_on + problem.params.line.built);

%% Make Sample
problem.samp = fraction_fact_samp(problem);

%% Output
problem.lock_on{z_idx} = lock_on;
problem.lock_off{z_idx} = lock_off;

end