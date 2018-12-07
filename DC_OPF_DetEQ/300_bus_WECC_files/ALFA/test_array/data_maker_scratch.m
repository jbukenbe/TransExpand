% This script was used to select the scenarios and their weights for
% importance sampleing, k means, and ALFA for use in the TNEP optimization.
% This scipt is not fully documented and does not automatically save the
% data to a file. The data sets that were saved for the TNEP paper are
% called is_data, km_data, lf_10_data, and lf_30_data.

% A large portion of this code is dedicated to checking the output
% approximations against a crossvalidation batch to ensure proper
% implimentation

%History            
%Version    Date        Who     Summary
%1          10/22/2018  JesseB  Documentation added *NOTE* originally created in June




%% Importance Sampeling
%{
scen_n_list = repelem([2 4 8 12 24],20);
scen_list = cell(100,1);
scen_w = cell(100,1);
global_min = min(scen_op_cost);
global_min = 0;
base = scen_op_cost - global_min;
base_m = mean(base);

p = 1/length(scen_op_cost);
q = (p*base)./(base_m);
cum_q = cumsum(q);
pi = p./q;


for s_idx = 1:100
    scens = zeros(1,scen_n_list(s_idx));

    for h_idx = 1:scen_n_list(s_idx)
        scens(h_idx) = find(rand() < cum_q,1);
    end
    
    scen_list{s_idx} = scens;
    scen_w{s_idx} = pi(scens)./scen_n_list(s_idx);
    
    
    
    apx_out = (output(:,scens)-global_min)*scen_w{s_idx}+global_min;
    [apx_out,t_id] = min(apx_out);
    mest(s_idx) = mean(output(t_id,:),2);
    
    
end

boxplot(mest,scen_n_list)

scen_n = repmat(repelem([1 5 10 20 30],1,5),1,3);
samp_lim = repelem([50 100 500],1,25);

scen_list = cell(75,1);
scen_w = cell(75,1);
global_mean = zeros(75,1);

imx = zeros(75,1);
for s_idx = 1:75
    [s_base,base_id] = sort(scen_op_cost,'descend');
    global_min = min(s_base);
    s_base = s_base - global_min;
    
    
    samp_n = scen_n(s_idx);
    samp_pool = samp_lim(s_idx);
    
    
    scen_list{s_idx} = base_id(randperm(samp_pool,samp_n));
    
    %scen_list{s_idx} = s_base(randperm(samp_pool,samp_n));
    scen_w{s_idx} = ones(1,samp_n)*samp_pool./(samp_n*8736);
    global_mean(s_idx) = mean(s_base(samp_pool:end))*(8736-samp_pool)/8736;

    imx(s_idx) = global_min + global_mean(s_idx) + scen_w{s_idx}*(scen_op_cost(scen_list{s_idx})-global_min);


end
plot(imx)
hold on



[sorted,sort_id] = sort(scen_op_cost,'descend');
ms = min(sorted);
sorted = sorted - ms;
adx = 500;
smpx = 20;
meen = mean(sorted(adx:end));
for s_idx = 1:300
    scen = output(s_idx,sort_id)-ms;
    
    imx(s_idx) = ms + (meen*(8736-adx)/8736) + sum(adx*scen(randperm(adx,smpx))./(smpx*8736));
    imp(s_idx) = ms +  mean(scen(randperm(8736,smpx)));
    
end
%plot(imx,'o')
%hold on
%plot(imp,'o')
tru = mean(output,2);

[simp,simp_id] = sort(imx);
plot(tru(simp_id))
hold on
plot(simp)
%}


%% K means Samples
%{
stoch_data = [pload.val;pRenew.val;pVarCost.val];
min_s = min(stoch_data,[],2);
stoch_u = stoch_data - min_s;
range_s = max(stoch_u,[],2);
range_s(range_s == 0) = 1;
stoch_u = stoch_u./range_s;
stoch_u = stoch_u';


scen_w = cell(100,1);
scen_load = cell(100,1);
scen_renew = cell(100,1);
scen_varcost = cell(100,1);

cluster_list = repelem([2 4 8 12 24],20);
for x_idx = 1:100
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



%}
%% Latent Factor Approximator Samples
run_n = 100;
plan_list = cell(run_n,1);
scen_list = cell(run_n,1);
global_mean = zeros(run_n,1);
scen_mean = cell(run_n,1);
scen_w = cell(run_n,1);
scen_w_pos = cell(run_n,1);

scen_n_list = repelem([2 4 8 12 24],20);
lf_n_list = repelem([2 3 7 10 15], 20);
plan_samp_size = 10;



% *NOTE* output is the hourly cost data that was output by the original
% sample of 300 networks, this was saved as 'output' but it is not the
% output of this script
[samp_set_size,~] = size(output);

true_cost = mean(output,2);
apx = zeros(samp_set_size,1);


x_idx = 1;
while x_idx < (run_n+1)
   %plan_samp_size = lf_plan_list(x_idx);
   scen_samp_size = scen_n_list(x_idx); 
   LF_n = lf_n_list(x_idx);
   
% pick plans
    plan_samp = randperm(samp_set_size, plan_samp_size);
    plan_op_costs = output(plan_samp,:);
    plan_list{x_idx} = plan_samp;

% problem mean
    global_mean(x_idx) = mean(mean(plan_op_costs));            
        
% pick scenarios
    A = plan_op_costs - mean(plan_op_costs);
    [U, S, V] = svd(A);
    best_scen = candexch(V(:,1:LF_n),scen_samp_size, 'tries',100,'display', 'off');
    
    % remove duplicates and add new scenarios up to required number
    bs = unique(best_scen);
    while ~isempty(best_scen)
        [cand,cand_id] = setdiff(V(:,1:LF_n),V(bs,1:LF_n),'rows');
        best_scen = candexch(cand,scen_samp_size-length(bs),'start',V(bs,1:LF_n), 'tries',100,'display', 'off');
        bs = unique([bs;cand_id(best_scen)]);
    end
    best_scen = bs;
    scen_list{x_idx} = best_scen;

% scenario means
    scen_mean{x_idx} = mean(plan_op_costs(:,best_scen));

% calculate raw scenario weights    
    X = (S(1:LF_n,1:LF_n)*V(best_scen,1:LF_n)')';
    H = V(:,1:LF_n)*S(1:LF_n,1:LF_n)*((X'*X)\(X'));    
    raw_w = sum(H);
    scen_w{x_idx} = raw_w;

% calculate positive weights with gradient descent
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

    else    
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
   
    % assign resulting positive weights
    if sum(sign(reg_w))~= length(reg_w)
        % process failed to produce positive weights, do not save this run
        scen_w_pos{x_idx}=-1;
    else
        scen_w_pos{x_idx} = reg_w;
        
        apx = ((output(:,best_scen) - scen_mean{x_idx})*scen_w_pos{x_idx}')/8736+global_mean(x_idx);
        sqer(x_idx) = (apx - true_cost)'*(apx - true_cost);
        
        x_idx = x_idx +1;
    end
    end
 
end



for t_idx = 1:100
        t_apx = ((output(:,scen_list{t_idx}) - scen_mean{t_idx})*scen_w_pos{t_idx}')/8736+global_mean(t_idx);
        t_sqer(t_idx) = (t_apx - true_cost)'*(t_apx - true_cost);  
        er_log(t_idx) = mean(t_apx-true_cost);
end

m = matfile('lf_10_data');
m.plan_list = plan_list;
m.scen_list = scen_list;
m.scen_mean = scen_mean;
m.scen_w = scen_w;
m.scen_w_pos = scen_w_pos;
m.global_mean = global_mean;
m.lf_n_list = lf_n_list;
m.scen_n_list = scen_n_list;
m.plan_samp_size = plan_samp_size;
m.sq_err = sqer;


%{
%% Test plotting

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
%}

    