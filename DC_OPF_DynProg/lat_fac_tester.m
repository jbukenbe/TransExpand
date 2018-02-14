

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

