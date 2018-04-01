function keep_lines = beta_graph_search(problem)
% This function takes the problem structure and a beta vector from a
% regression to use graph search methods to return lines that appear in low
% cost and reliable plans

%History            
%Version    Date        Who     Summary
%1          12/02/2017  JesseB  Initial Version adapted from beta_explorer

%% Initialize data
% get problem data
beta = problem.pls.beta;
good_plan_threshold = problem.pls.good_plan_threshold;
line_id = problem.params.cand.line_id;
line_n = length(line_id);
line_cost = problem.params.new_line_cost;
line_beta = beta(2:size(line_id,1)+1);

% make mapping from beta coefficients to lines;
interact_list = nchoosek(line_id,2);
beta_map = [0 0;line_id, zeros(size(line_id,1),1);interact_list];
beta_line_map = [beta, beta_map];


for l_idx = 1:length(line_id)
    line_idx = line_id(l_idx);
    top_b_logic = beta_line_map(:,2) == line_idx;
    bot_b_logic = beta_line_map(:,3) == line_idx;
    line_beta(l_idx) = beta_line_map(:,1)'*(top_b_logic+bot_b_logic);
end


%% Run A* Graph Search
% set starting position
this_plan = zeros(1,line_n);

% set list of visited states and possible transitions
visited_states = bi2de(this_plan);
new_transitions = eye(line_n);
new_trans_n = size(new_transitions,1);
new_transitions_id = bi2de(new_transitions);
transition_list = [];
good_plan_id = [];
good_plan_cost = [];
expected_costs = [];

% loop through list of possible transitions looking for good plans
stop = 0; i_idx = 1;
while ~stop
    % find transition costs and fits
    transition_costs = new_transitions*line_cost;
    transition_fits = x2fx(new_transitions,'interactions')*beta;
    
    % find any plans with acceptable expected operating costs
    new_good_plan_logic = (transition_fits < good_plan_threshold);
    good_plan_id = [good_plan_id; new_transitions_id(new_good_plan_logic)];
    good_plan_cost = [good_plan_cost; transition_costs(new_good_plan_logic)];
    
    % remove good plans
    transition_fits(new_good_plan_logic) = [];
    new_transitions(new_good_plan_logic,:) = [];
    transition_costs(new_good_plan_logic) = [];
    new_trans_n = new_trans_n - sum(new_good_plan_logic);
    if length(good_plan_id) > 10000 || (i_idx > 5000)
        stop = 1;
    else
        i_idx = i_idx +1;
        % update beta for next line heuristic
        beta_h = repmat(line_beta',new_trans_n,1);
    
        % find expected line costs for other transition plans
        beta_rate = beta_h./line_cost';
        improve_rate = max(-beta_rate,[],2);
        over_op_cost = 20*(transition_fits - good_plan_threshold);
        transition_h = over_op_cost./improve_rate;
        transition_e_cost = transition_costs + transition_h;
        
        % update global problem lists;
        expected_costs = [expected_costs; transition_e_cost];
        transition_list = [transition_list; new_transitions_id];
        
        % pick new state based on minimum expected final cost
        [~, new_state_index] = min(expected_costs);
        visited_states = [visited_states; transition_list(new_state_index)];
        this_plan = de2bi(transition_list(new_state_index),line_n);
        expected_costs(new_state_index) = [];
        transition_list(new_state_index) = [];
        
        % find new transition states
        new_transitions = eye(line_n) + this_plan;
        new_transitions(max(new_transitions,[],2)>1,:)=[];
        
        %remove known transition states
        new_transitions_id = bi2de(new_transitions);
        new_transitions_id = setdiff(new_transitions_id,[transition_list;visited_states]);
        new_transitions = de2bi(new_transitions_id,line_n);
        new_trans_n = size(new_transitions,1);
    end
end

% find lines to keep and those to drop
[~,best_line_id] = sort(sum(de2bi(good_plan_id,line_n)),'descend');
keep_lines = line_id(best_line_id(1:ceil(line_n/2)));

end