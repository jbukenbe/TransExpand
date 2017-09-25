function problem = pca_val_est(problem)
% This function takes the problem from the dyn_DC_OPF function and uses
% principal component analysis to make estimates of the value of future
% plans 

r_log = zeros(1,300);
for t_idx = 1:300
all_plan = de2bi(problem.plan_id-1, problem.params.cand.n);
ordering = randperm(1024)';
plan = all_plan(ordering(1:500),:);
val = problem.cand_cost(ordering(1:500));

A = [val,plan];
Amean = mean(A);
Astd = std(A);
A_stand = (A-Amean)./Astd;

[eiganvectors,~, eiganvalues] = pca(A_stand);
for e_idx = 1: size(eiganvalues)-1
    line_trans = all_plan*eiganvectors(2:end,1:e_idx);
    stand_val_est = (line_trans)*eiganvectors(1,1:e_idx)';
    val_est = stand_val_est*Astd(1)+Amean(1);
    mdl = fitlm(val_est,problem.cand_cost);
    if e_idx == 8
        r_log(t_idx) = mdl.Rsquared.Ordinary;
    end
    %scatter(val_est, problem.cand_cost);
    %xlabel('PCA Estimate');
    %ylabel('True Value')
end
end
    r_log'
    me = mean(r_log,2);
    mx = max(r_log')';
    mn = min(r_log')';
    s=std(r_log')';
    pca_data.val_est = val_est;
    pca_data.eiganvalues = eiganvalues;
    pca_data.eiganvectors = eiganvectors;
    problem.params.pca_data =  pca_data;



end


