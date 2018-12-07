classdef CheckTSTNEPClass
% This class is used to check the network solution to the tstnep problem 
        
%History            
%Version    Date        Who     Summary
%1          11/30/2018  JesseB  Initial version
%
%
%   TODO:   Plot results method
%           Better handling of generation plan data
%           
%
    
    properties       
        cluster_n
        cluster_w
        cluster_mean
        cluster_list
        
        apx_samp_n       
        apx_mean
        apx_w
        
        run_n
        hr_list
        gn_list
        gen_in
        
        texp_plan
        line_in
        opf_results
        
        cluster_op_cost 
        cluster_line_cost
        total_clust_cost
        total_weighted_cost
        
    end
    
    methods
        function obj = CheckTSTNEPClass(Luca, gexp_full, apx_samp_n)
            % initialization function
            if nargin < 3
                obj.apx_samp_n = 5;
            else
                obj.apx_samp_n = apx_samp_n;
            end
                        
            % take cluster set from LUCA for evaluation
            obj.cluster_n = Luca.cluster_n;

            % get new sample gen hours for each cluster
            obj.run_n = obj.apx_samp_n*obj.cluster_n;
            [~, Luca_lexp_samp_out] = Luca.get_lexp_samp(obj.apx_samp_n);

            % store sample data
            obj.cluster_list = Luca_lexp_samp_out.lexp_clust_samp;
            obj.cluster_mean = Luca_lexp_samp_out.cluster_mean;
            obj.cluster_w = Luca_lexp_samp_out.cluster_w;
            obj.apx_mean = Luca_lexp_samp_out.lexp_samp_mean;
            obj.apx_w = Luca_lexp_samp_out.lexp_w;
            
            obj.hr_list = Luca_lexp_samp_out.lexp_hr_samp';
            obj.gn_list = Luca_lexp_samp_out.lexp_gexp_samp';
            
            % format hr and gen list
            obj.hr_list = obj.hr_list(:);

            obj.gn_list = obj.gn_list(:);
            obj.gen_in = gexp_full(obj.gn_list,:);
% TODO: only store unique values and reference properly later
        end
        
        
        function obj = get_opf_runs(obj, texp_plan)
            % take an input transmission plan and calculate the OPFs needed
            % for a solution
            
            % store plan
            obj.texp_plan = texp_plan;
            
            % get lines used in each cluster
            obj.line_in = zeros(obj.cluster_n,51);
            line_out = [texp_plan(:,4), texp_plan(:,1)];
            lin_idx = sub2ind(size(obj.line_in),line_out(:,1), line_out(:, 2));
            obj.line_in(lin_idx) = 1;

            % make opf object 
            solution_check = OptOPFClass(obj.hr_list, obj.line_in, obj.gen_in);

            % run opf problems
            output = zeros(obj.run_n,1);
            for r_idx = 1:obj.run_n
                hr_idx = r_idx;
                ln_idx = obj.cluster_list(r_idx);
                gn_idx = r_idx;
                
                [~,output(r_idx)] = solution_check.solve(hr_idx, ln_idx, gn_idx);
            end

            obj.opf_results = output;
    
        end
        
        
        function obj = check_plan(obj, line_cost)           
            % approximate costs of transmission network expansion plan
            opf_out_table = reshape(obj.opf_results, obj.apx_samp_n, obj.cluster_n)';
            obj.cluster_op_cost = sum((opf_out_table - obj.apx_mean).*obj.apx_w,2) + obj.cluster_mean;
            obj.cluster_line_cost = obj.line_in*line_cost;

            obj.total_clust_cost = obj.cluster_op_cost + obj.cluster_line_cost;
            obj.total_weighted_cost = obj.total_clust_cost'*obj.cluster_w;   
        end
        
        
        function obj = plot_results(obj)
            % plots the cost results of the tstnep plan
            
            
        end
    end
end


