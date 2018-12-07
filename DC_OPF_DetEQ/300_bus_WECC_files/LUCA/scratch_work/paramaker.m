param_set = cell(8,1);
load tstnep_run_3 params
p_in = params;
clear params
for s_idx = 2:9
    idx = (1+4*(s_idx-1)):(s_idx)*4;
    
    sn_rhm = [p_in.hr_run_map(1:4); p_in.hr_run_map(idx)];
    sn_grm = [p_in.gexp_run_map(1:4); p_in.gexp_run_map(idx)];
    sn_crm = [p_in.clust_run_map(1:4); 2*ones(4,1)];
    sn_cn = 2;
    sn_bc = p_in.baseline_clust;

    sn_gl = [p_in.gexp_gens(1:4,:); p_in.gexp_gens(idx,:)];
    sn_cw = [.5;.5];
    sn_wp = [p_in.run_w_pos(1:4); p_in.run_w_pos(idx)];
    sn_w = [p_in.run_w(1:4); p_in.run_w(idx)];
    sn_rmc = [p_in.run_mean_center(1:4); p_in.run_mean_center(idx)];
    
    sn_cm = [p_in.cluster_mean(1); p_in.cluster_mean(s_idx)];
    sn_cs = p_in.gexp_cluster_save;
    
    
 
    
    params.hr_run_map = sn_rhm;
    params.gexp_run_map = sn_grm;        
    params.clust_run_map = sn_crm;
    params.cluster_n = sn_cn;
    params.baseline_clust = sn_bc;
                   
    params.gexp_gens = sn_gl;
    params.cluster_w = sn_cw;
    params.run_w_pos = sn_wp;
    params.run_w = sn_w;
    params.run_mean_center = sn_rmc;            
    params.cluster_mean = sn_cm;
    params.gexp_cluster_save = sn_cs; 
  
  
    param_set{s_idx-1} = params;
end



m = matfile('tstnep_param_data2');
m.param_set = param_set;