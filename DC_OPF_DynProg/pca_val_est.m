function problem = pca_val_est(problem)
% This function takes the problem from the dyn_DC_OPF function and uses
% principal component analysis to make estimates of the value of future
% plans 


val = problem.cand_cost;
plan = de2bi(problem.plan_id-1, problem.params.cand.n);
A = [val,plan];
Amean = mean(A);
Astd = std(A);
A_stand = (A-Amean)./Astd;

[eiganvectors,~, eiganvalues] = pca(A_stand);

line_trans = plan*eiganvectors(2:end,:);
stand_val_est = (line_trans.*eiganvalues')*eiganvectors(1,:)';
val_est = stand_val_est*Astd(1)+Amean(1);


pca_data.val_est = val_est;
pca_data.eiganvalues = eiganvalues;
pca_data.eiganvectors = eiganvectors;
problem.params.pca_data =  pca_data;



end


