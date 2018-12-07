function [] = lf_w_tuner(output,base,cross_val)

run_n = 100;
plan_list = cell(run_n,1);
scen_list = cell(run_n,1);
global_mean = zeros(run_n,1);
scen_mean = cell(run_n,1);
scen_w = cell(run_n,1);
scen_w_pos = cell(run_n,1);

%scen_n_list = repelem([2 4 8 12 24],20);
%lf_n_list = repelem([2 4 7 13 13], 20);


scen_n_list = repelem([2 4 8 12 24],20);
lf_n_list = repelem([2 4 7 13 13], 20);
plan_samp_size = 10;

true_cost = mean(output,2);
cross_cost = mean(cross_val,2);
[~,cross_sort_id] = sort(cross_cost);

x_idx = 1;
retry = 20;
fail = 0;
while x_idx < (run_n+1)
   %plan_samp_size = lf_plan_list(x_idx);
   scen_samp_size = scen_n_list(x_idx); 
   LF_n = lf_n_list(x_idx);
   
% pick plans
    plan_samp = zeros(plan_samp_size,1);
    [~,sort_id] =sort(mean(output,2));
    for d_idx = 1:plan_samp_size
        range_low = 1+(d_idx-1)*(300/plan_samp_size);
        range_hi = d_idx*(300/plan_samp_size);
        plan_samp(d_idx) = sort_id(randi([range_low, range_hi]));
    end
    %plan_samp = randperm(300, plan_samp_size);
    plan_op_costs = output(plan_samp,:);
    plan_list{x_idx} = plan_samp;

% problem mean
    global_mean(x_idx) = mean(mean(plan_op_costs));            
        
% pick scenarios
    A = plan_op_costs - mean(plan_op_costs);
    [~, S, V] = svd(A);
    best_scen = candexch(V(:,1:LF_n),scen_samp_size, 'tries',100,'display', 'off');
    bs = unique(best_scen);
 %   while ~isempty(best_scen)
 %       [cand,cand_id] = setdiff(V(:,1:LF_n),V(bs,1:LF_n),'rows');
 %       best_scen = candexch(cand,scen_samp_size-length(bs),'start',V(bs,1:LF_n), 'tries',100,'display', 'off');
 %       bs = unique([bs;cand_id(best_scen)]);
 %   end
 %   best_scen = bs;

    
base_m = mean(base);

p = 1/length(base);
q = (p*base)./(base_m);
cum_q = cumsum(q);
pi = p./q; 
    

    for h_idx = 1:scen_n_list(x_idx)
        is_scens(h_idx) = find(rand() < cum_q,1);
    end    
    
is_w =  pi(is_scens)./scen_n_list(x_idx);   

scens_needed = scen_n_list(x_idx)-length(bs);
rand_scen = is_scens(randperm(scens_needed));
best_scen = [bs;rand_scen'];

lfis_w = pi(best_scen)./length(best_scen);    
lfis_scens = best_scen;

is_big_apx = output(:,is_scens)*is_w;
is_fine_apx = cross_val(:,is_scens)*is_w;
lfis_big_apx = output(:,lfis_scens)*lfis_w;
lfis_fine_apx = cross_val(:,lfis_scens)*lfis_w;

hold off
plot(true_cost(sort_id))
hold on
plot(cross_cost(cross_sort_id))
plot(is_big_apx(sort_id))
plot(is_fine_apx(cross_sort_id))
plot(lfis_big_apx(sort_id))
plot(lfis_fine_apx(cross_sort_id))



    
    
    
    
 %{   
    scen_list{x_idx} = best_scen;

% scenario means
    scen_mean{x_idx} = mean(plan_op_costs(:,best_scen));

% calculate raw scenario weights    
    X = (S(1:LF_n,1:LF_n)*V(best_scen,1:LF_n)')';
    H = V(:,1:LF_n)*S(1:LF_n,1:LF_n)*((X'*X)\(X'));
    
    H_real = H;
    H_real(best_scen,:) = eye(scen_samp_size);
    raw_w_real = sum(H_real);
    
    raw_w = sum(H);
    scen_w{x_idx} = raw_w;


 %calculate positive weights with penalized gradient descent
    reg_w = raw_w;
    neg_w = find(sign(raw_w) < 0);
    pos_w = find(sign(raw_w) > 0);
    target_cost = sum(A,2);
    A_sub = A(:,best_scen);
    count = 0;

    if isempty(neg_w)
        scen_w_pos{x_idx} = reg_w;
        
        apx = ((output(:,best_scen) - scen_mean{x_idx})*scen_w_pos{x_idx}')/8736+global_mean(x_idx);
        sqer(x_idx) = (apx - true_cost)'*(apx - true_cost);
        
        x_idx = x_idx +1;
        retry = 0;
    elseif retry < max(12 - LF_n,0)
        retry = retry + 1
    else
        
 %{   
    fail = fail + 1
    % continue until no negative weights exist or convergence error 
    while (~isempty(neg_w)) && (count < 750)
        %penalize negative weights
        reg_w(neg_w) = reg_w(neg_w) +abs(raw_w(neg_w)*.02*(count/100));
        
        % adjust positive weights to accommodate inflating negative values
        for d_idx = 1:100
            err = A_sub*reg_w' - target_cost;
            row_order = randperm(plan_samp_size);
            alpha = max(abs(reg_w))/(max(count,100)*max(err)*max(sum(abs(A_sub(:,pos_w)),2)));
            for r_idx = row_order
                reg_w(pos_w) = reg_w(pos_w) - alpha*err(r_idx)*A_sub(r_idx,pos_w);
            end
        end
        
        % check for positive and negative weights
        neg_w = find(sign(reg_w) < 0);
        pos_w = find(sign(reg_w) > 0);
        count = count+1;
    end
   
 %}   
    
 
    
    
    
    
    
    
    reg_w = raw_w;
    add_w = raw_w*0;
    reg_gm = global_mean(x_idx);
    
    neg_w = find(sign(raw_w) < 0);
    pos_w = find(sign(raw_w) > 0);
    target_cost = mean(plan_op_costs,2);
    A_sub = A(:,best_scen);
    Z_sub = plan_op_costs(:,best_scen);
    
    apx_fix = A_sub*raw_w';
    
    
    count = 0;
    
        % continue until no negative weights exist or convergence error 
    while (~isempty(neg_w)) && (count < 750)
        % Deflate global mean
        reg_gm = reg_gm*.998;
        
        %penalize negative weights
        %reg_w(neg_w) = reg_w(neg_w) +abs(raw_w(neg_w)*.02*(count/100));
        
        % adjust positive weights to accommodate inflating negative values
        for d_idx = 1:300
            err = (A_sub*reg_w' + Z_sub*add_w')./8736 + reg_gm - target_cost;
            row_order = randperm(plan_samp_size-10);
            
            alpha = min(abs(raw_w))/((count+100)*max(abs(err))*max(sum(abs(Z_sub),2)));
            %omega = max(abs(reg_w))/(max(count,100)*max(err)*max(sum(abs(A_sub(:,pos_w)),2)));
            for r_idx = row_order
                add_w = add_w - alpha*err(r_idx)*Z_sub(r_idx,:);
             %   reg_w(pos_w) = reg_w(pos_w) - omega*err(r_idx)*A_sub(r_idx,pos_w);
            end
        end
        
        % check for positive and negative weights
        neg_w = find(sign(reg_w + add_w) < 0);
        pos_w = find(sign(reg_w + add_w) > 0);
        count = count+1;
    end
    
    
    
    
    
  
    
    
    
    % assign resulting positive weights
    if sum(sign(reg_w))~= length(reg_w)
        scen_w_pos{x_idx}=-1
    else
        scen_w_pos{x_idx} = reg_w;
        
        apx = ((output(:,best_scen) - scen_mean{x_idx})*scen_w_pos{x_idx}')/8736+global_mean(x_idx);
        sqer(x_idx) = (apx - true_cost)'*(apx - true_cost);
        
        x_idx = x_idx +1;
        retry = 0;
    end
    end
 %}
end



for t_idx = 1:100
        t_apx = ((output(:,scen_list{t_idx}) - scen_mean{t_idx})*scen_w_pos{t_idx}')/8736+global_mean(t_idx);
        t_sqer(t_idx) = (t_apx - true_cost)'*(t_apx - true_cost);  
        er_log(t_idx) = mean(t_apx-true_cost);
end





end