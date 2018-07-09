%% K means Samples
%stoch_data = [pload.val;pRenew.val;pVarCost.val];
%min_s = min(stoch_data,[],2);
%stoch_u = stoch_data - min_s;
%range_s = max(stoch_u,[],2);
%range_s(range_s == 0) = 1;
%stoch_u = stoch_u./range_s;
%stoch_u = stoch_u';

%{
scen_w = cell(60,1);
scen_load = cell(60,1);
scen_renew = cell(60,1);
scen_varcost = cell(60,1);

cluster_list = repelem(1:30,2);
for x_idx = 1:60
    cluster_n = cluster_list(x_idx);
    [idx, C] = kmeans(stoch_u,cluster_n);
    
    C = C'.*range_s + min_s; 
    
    scen_load{x_idx} = C(1:312,:);
    scen_renew{x_idx} = C(313:961,:);
    scen_varcost{x_idx} = C(962:1941,:);
    
    w = zeros(cluster_n,1);
    for c_idx = 1:cluster_n
        w(c_idx) = nnz(idx == c_idx);
    end
    scen_w{x_idx} = w;
    clear w
end

%% Monte Carlo Samples

scen_list = cell(25,1);
scen_w = cell(25,1);

mc_size_list = repelem([1;5;10;20;30],5);
for x_idx = 1:length(mc_size_list)
   samp_size = mc_size_list(x_idx); 
   scen_list{x_idx} = randperm(8736, samp_size);
   scen_w{x_idx} = ones(samp_size,1)*8760/samp_size;
end

%% Latent Factor Approximator Samples
%plan_list = cell(75,1);
%scen_list = cell(75,1);
global_mean = zeros(75,1);
scen_mean = cell(75,1);
scen_w = cell(75,1);
scen_w_pos = cell(75,1);

lf_scen_size_list = [1; 5; 10; 20; 30];
lf_scen_size_list = repelem(lf_scen_size_list, [20,20,15,10,10]);
lf_plan_count = [5; 10; 30; 50];
lf_plan_list = [repelem(lf_plan_count,5);
                repelem(lf_plan_count,5);
                repelem(lf_plan_count(2:4),5);
                repelem(lf_plan_count(3:4),5);
                repelem(lf_plan_count(3:4),5)];

redo = [43,68];
for x_idx = redo
   plan_samp_size = lf_plan_list(x_idx);
   scen_samp_size = lf_scen_size_list(x_idx); 
   
% pick plans
    plan_samp = randperm(300, plan_samp_size);
    plan_op_costs = output(plan_samp,:);
    plan_list{x_idx} = plan_samp;

% problem mean
    global_mean(x_idx) = mean(mean(plan_op_costs));            
        
% pick scenarios
    A = plan_op_costs - mean(plan_op_costs);
    [~, S, V] = svd(A);
    LF_n = max(1,scen_samp_size-1);
    best_scen = candexch(V(:,1:LF_n),scen_samp_size, 'tries',100,'display', 'off');
    scen_list{x_idx} = best_scen;

% scenario means
    scen_mean{x_idx} = mean(plan_op_costs(:,best_scen));

% calculate raw scenario weights    
    X = (S(1:LF_n,1:LF_n)*V(best_scen,1:LF_n)')';
    H = V(:,1:LF_n)*S(1:LF_n,1:LF_n)*((X'*X)\(X'));    
    raw_w = sum(H);
    scen_w{x_idx,1} = raw_w;

% calculate positive weights with penalized gradient descent
    reg_w = raw_w;
    neg_w = find(sign(raw_w) < 0);
    pos_w = find(sign(raw_w) > 0);
    target_cost = sum(A,2);
    A_sub = A(:,best_scen);
    count = 0;

    % continue until no negative weights exist or convergence error 
    while (~isempty(neg_w)) && (count < 750)
        %penalize negative weights
        reg_w(neg_w) = reg_w(neg_w) +abs(raw_w(neg_w)*.02*(count/100));
        
        % adjust positive weights to accommodate inflating negative values
        for d_idx = 1:100
            err = A_sub*reg_w' - target_cost;
            row_order = randperm(scen_samp_size);
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
    
    % assign resulting positive weights
    if sum(sign(reg_w))~= length(reg_w)
        scen_w_pos{x_idx}=-1
    else
        scen_w_pos{x_idx} = reg_w;
    end
    

end



%% Test plotting
%}
actual = mean(output,2);
[sor,id]= sort(actual);
for x_idx = 1:75
    plan_id = plan_list{x_idx};
    scens = scen_list_lf{x_idx};
    scen_m = scen_mean{x_idx};
    gmean = global_mean(x_idx);
    sw = scen_w_lf{x_idx};
    swp = scen_w_pos{x_idx};
    
    swa = ((output(:,scens)-scen_m)*sw')./8760+gmean;
    swpa = ((output(:,scens)-scen_m)*swp')./8760+gmean;
    
    mc_idx = min(x_idx,25);
    scens = scen_list{mc_idx};
    mcw = scen_w{mc_idx};
    mca = (output(:,scens)*mcw)./8760;
    
    
    km_idx = min(x_idx,60);
    c_n = length(km_w{km_idx});
    c_vec = [scen_load{km_idx};scen_renew{km_idx};scen_varcost{km_idx}];
    sckm = zeros(c_n,1);
    for c_idx = 1:c_n
        di = zeros(8736,1);
        for h_idx = 1:8736
            di(h_idx) = norm(c_vec(:,c_idx)-stoch_data(:,h_idx),2);
        end
        [~,sckm(c_idx)] = min(di);
    end
    kma = (output(:,sckm)*km_w{km_idx})./8760;
    
    plot(sor)
    hold on
    plot(swa(id));
    plot(swpa(id));
    plot(mca(id));
    plot(kma(id));
    hold off
legend('real','raw','positive','mc','km');
end


    