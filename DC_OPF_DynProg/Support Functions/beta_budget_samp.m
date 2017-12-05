function samp = beta_budget_samp(problem)
% This function returns a sample based on the line budjet and beta
% information for current lines

%History            
%Version    Date        Who     Summary
%1          12/4/2017  JesseB  Initial Version

% get required problem data
beta = problem.pls.beta{problem.z_idx-1};
cand_n = problem.params.cand.n(problem.z_idx-1);
samp_n = problem.params.refine_samp_n;
line_costs = problem.params.new_line_cost{problem.z_idx-1};
best_budget = (min(problem.cand_full_cost)- min(problem.cand_op_cost));
budget_limit_hi = 1.2*best_budget;
budget_limit_low = (1-.618)*best_budget;
plan_threshold = min(problem.cand_op_cost)*1.3;

% get lines to highlight this search
[~,beta_order] = sort(beta(2:end));
keep_lines = beta_order(1:10);


% initialize sample storage
samp = zeros(samp_n,cand_n,'uint8');
samp_idx = 1;

% find sample
while samp_idx <= samp_n
    %initialize loop data
    this_plan = zeros(1,cand_n); 
    plan_budget = 0;
    plan_fit = 0;
    l_idx = 1;
    stop = 0;
    
    % get first two highlighted lines
    good_line = randperm(10,1);
    this_plan(keep_lines(good_line)) = 1;
    
    % randomly generate potential ordering of lines
    line_order = randperm(cand_n);

    while ~stop 
        % add one line to this plan
        this_plan(line_order(l_idx)) = 1;
        plan_budget = plan_budget + line_costs(line_order(l_idx));
        
        % estimate plan fit
        if problem.params.pls.interaction
            this_plan_interaction = x2fx(this_plan,'interactions');
            plan_fit = this_plan_interaction*beta;
        else
            plan_fit = [1,this_plan]*beta; 
        end
        
        
        if plan_budget < budget_limit_low
        % keep adding lines until a reasonable budget is reached
            l_idx = l_idx +1;
        elseif plan_budget > budget_limit_hi
        % if line budget is exceeded start with new line ordering
            stop = 1;
        elseif plan_fit < plan_threshold
        % if the plan is not over budget but should perform well add to sample
            samp(samp_idx,:) = uint8(this_plan);
            samp_idx = samp_idx + 1;
            stop = 1;
        else
        % otherwise add another line 
            l_idx = l_idx + 1;
        end
    end

end
end
    
   