function problem = pca_val_est(problem)
% This function takes the problem from the dyn_DC_OPF function and uses
% principal component analysis to make estimates of the value of future
% plans 

%History            
%Version    Date        Who     Summary
%1          09/17/2017  JesseB  Initial Version

plan_log = zeros(30,10);
all_plan = de2bi(problem.plan_id-1, problem.params.cand.n);
ordering = randperm(1024)';
plan = all_plan(ordering(1:30),:);
val = problem.cand_cost(ordering(1:30));

for t_idx = 1:10

A = [val,plan];
Amean = mean(A);
Astd = std(A)+eps;
A_stand = (A-Amean)./Astd;

[eiganvectors,~, eiganvalues] = pca(A_stand);
line_trans = all_plan*eiganvectors(2:end,1:end-3);
stand_val_est = (line_trans)*eiganvectors(1,1:end-3)';
val_est = stand_val_est*Astd(1)+Amean(1);

mdl = fitlm(val_est,problem.cand_cost);
scatter(val_est, problem.cand_cost);
xlabel('PCA Estimate');
ylabel('True Value')

[~, plan_id] = sort(val_est);
plan_log(:,t_idx) = problem.cand_cost(plan_id(1:30));
val = problem.cand_cost(plan_id(1:30));
plan = all_plan(plan_id(1:30),:);

end
    me(t_idx) = mean(plan_log);
    mx(t_idx) = max(plan_log);
    mn(t_ids) = min(plan_log);

    pca_data.val_est = val_est;
    pca_data.eiganvalues = eiganvalues;
    pca_data.eiganvectors = eiganvectors;
    problem.params.pca_data =  pca_data;



end


