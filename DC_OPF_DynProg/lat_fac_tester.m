
T = [A(1:20,:);A(32238:32273,:);A(randperm(32768,100),:)];
Tavg = mean(T);
T_mean_cent = T - Tavg;
[U,S, V] = svd(T_mean_cent);
X = (S(1:3,1:3)*[V(2:3,1:3);V(8,1:3)]')';
H = V(:,1:3)*S(1:3,1:3)*inv(X'*X)*X';
Hridge = V(:,1:3)*inv(X'*X -.001*eye(3))*X';
sumH = sum(H);
sumH_ridge = sum(Hridge);
for a_idx = 1:600
    partial_vec = A(a_idx+22354,:)- Tavg;
    partial_vec(1) = 0;
    partial_vec(4:7) = 0;
    partial_vec(9) = 0;
    Y = partial_vec(logical([0 1 1 0 0 0 0 1 0]))';
    [approx_out, u_approx] = mat_fill_svd(S(1:3,1:3), V,partial_vec);
    approx_out = approx_out+ Tavg-A(a_idx+300,:)
    y_hat = (H*Y)'+ Tavg - A(a_idx+22354,:)
    y_hat_ridge = (Hridge*Y)'+ Tavg - A(a_idx+22354,:)
    expect = sumH*Y/9 + mean(Tavg)
    expect_ridge = sumH_ridge*Y/9 + mean(Tavg)
    tr = mean(A(a_idx+22354,:))
end




%{
robust_log = zeros(20,20);
expected_cost_err_log = zeros(20,20);

for svd_scen = 1:20
    for init_idx = 10:10:200
        test_prob = dyn_pls_demo(saved_problem, init_idx+25,svd_scen);
        second_samp_approx = test_prob.scen_op_cost((26+init_idx):end,:);
        
        robust_error = (abs(max(saved_problem.second_samp_real(36:end,:)')-max(second_samp_approx')))./max(saved_problem.second_samp_real(36:end,:)');
        robust_mse = sum(robust_error)/3;
        robust_log(init_idx/10,svd_scen) = robust_mse;
                
        
        second_samp_err = (second_samp_approx - saved_problem.second_samp_real(36:end,:));
        second_samp_sq_err = second_samp_err.*second_samp_err;
        second_samp_mse = sum(sum(second_samp_sq_err))/(300*103);
        expected_cost_err_log(init_idx/10,svd_scen) = second_samp_mse;
    end
end

%}