classdef LucaClass
% This is the Latent Uncertainty Clustering Algorithm (LUCA) class. The
% class lets an instance of LUCA be created quickly so different data sets
% and parameters can easily be compared
    
    
%History            
%Version    Date        Who     Summary
%1          11/13/2018  JesseB  Created (not working)
%2          11/19/2018  JesseB  HOSVD and gexp sample working
%3          11/23/2018  JesseB  Fixed bug with gexp samp index
%4          11/27/2018  JesseB  Full working version (added Uhat, cluster, 
%                               TSTNEP samp, and gradient decent methods)
%5          11/28/2018  JesseB  Improved TSTNEP cluster with cluster_hr_means
%                               Added: Make TSTNEP parameters
%                               Changed mean centering behavior
%6          11/29/2018  JesseB  Moved first scenario to own cluster
%
%
%
% TODO:     Better use of tol or dim
%           D and G efficiency measurements
%           Cluster efficiency measurements
%           


    properties
        Z
        size_z
        Ugn
        Sgn
        Vgn

        params
        
        gexp_samp_id
        gexp_samp
        gexp_samp_n
        gexp_X
        gexp_mean
        Uhat_gexp
        
        cluster_n
        cluster_cen
        gexp_cluster
        cluster_w
        cluster_sz
        
        cluster_mean
        lexp_samp_n
        lexp_samp_id
        lexp_hr_samp
        lexp_gexp_samp
        lexp_clust_samp
        lexp_samp_mean
        lexp_w
        lexp_wpos
    end
    
    methods
        function obj = LucaClass(data_in, params)
            % Defines the Luca object and saves the mean centered data
            obj.Z = data_in;
            obj.size_z = size(data_in);
            
%[z_unfold, obj.size_z] = obj.unfold(1, data_in);
%obj.hr_means = mean(z_unfold,2);

%obj.A = obj.fold(1, obj.size_z, (z_unfold - obj.hr_means));
            
            % Save LUCA params if none are specified for each param
            default_params = struct(    'tol', .01,...
                                        'dim', obj.size_z-1,...
                                        'verbose', 0);                                
            if nargin < 2
                params = {};
            end                                    
            if isempty(params)
                obj.params = default_params;
            else
                fields = fieldnames(default_params);
                for i = 1:length(fields)
                    if ~isfield(params,fields{i})
                        params.(fields{i}) = default_params.(fields{i});
                    end
                end
                obj.params = params;
            end
        end    
        
        
        function [obj, gexp_samp, gexp_X, best_samp] = get_gexp_samp(obj, samp_n)
            % Generate a sample with size samp_n of (hour, line)
            % coordinates and the design matrix X that will be used to
            % calculate the gen beta values for unobserved generation plans
            
            % record desired sample size
            obj.gexp_samp_n = samp_n;
            
            
            % latent factors that can be approximated
            if samp_n <= obj.params.dim(2)
                lf_n = samp_n - 1;
            else
                lf_n = obj.params.dim(2);
            end
            
            
            % get SVD info            
            Zgn = obj.unfold(2, obj.Z);
            Zgn_mean = mean(Zgn);
            [obj.Ugn, obj.Sgn, obj.Vgn] = svd(Zgn-Zgn_mean,'econ');

            
            % retrieve scaled right singular vectors           
            Xgn = (obj.Sgn*(obj.Vgn'))';
            Xgn = Xgn(:,1:lf_n);
           
            
            % pick near D-optimal line_hr sample to estimate Ugn_hat
            [gexp_X, best_samp]= obj.D_opt_samp(Xgn, samp_n, 2);
            obj.gexp_X = gexp_X;
            obj.gexp_mean = Zgn_mean(best_samp);                

            
            % convert best_samp to lexp and hour IDs
            gexp_samp_tl = mod(best_samp-1, obj.size_z(3))+1;
            gexp_samp_hr = ceil(best_samp./obj.size_z(3));
            gexp_samp = [gexp_samp_hr, gexp_samp_tl];
            
            % output sample
            obj.gexp_samp_id = best_samp;
            obj.gexp_samp = gexp_samp;
        end
        
        
        function obj = get_uhat_gexp(obj, samp_in)
            % Estimates the value of U for gexp data based on input sample
            
            % reshape data into table with hr and gens for mean centering
            samp_in_n = numel(samp_in);
            samp_in_table = reshape(samp_in, samp_in_n/obj.gexp_samp_n, obj.gexp_samp_n);
            
            % mean center sample data
            mc_samp = samp_in_table - obj.gexp_mean;
            
            % approximate beta for generation plans
            beta_maker = (obj.gexp_X'*obj.gexp_X)\(obj.gexp_X');            
            obj.Uhat_gexp = (beta_maker*mc_samp')';
            
            % create debug plots
            if obj.params.verbose
                obj.plot_Uhat_gexp(beta_maker)
            end
        end
        
        
        function obj = cluster_gexp(obj, cluster_n, lf_n)
            % Group the generation plans into cluster_n clusters based on
            % the similarity of the first lf_n vectors in Uhat_gexp 
            
            
            if nargin < 3
                lf_n = obj.params.dim(2);
                if nargin < 2
                    cluster_n = lf_n;
                end
            end 
            obj.cluster_n = cluster_n;
            
            % scale Uhat data to proportion
            Uhat_clust = obj.Uhat_gexp(:,1:lf_n)*obj.Sgn(1:lf_n, 1:lf_n);
            
            % make clusters removing the baseline scenario
            [clust_ex_baseline, clust_cen_ex_baseline] = kmeans(Uhat_clust(2:end,:), cluster_n-1, 'dist','sqeuclidean', 'replicates',10);
            
            % keep the baseline as cluster 1
            obj.gexp_cluster = [1; clust_ex_baseline + 1];
            obj.cluster_cen = [Uhat_clust(1,:);clust_cen_ex_baseline];
            
            % save cluster weights
            obj.cluster_sz = zeros(cluster_n,1);
            obj.cluster_w = zeros(cluster_n,1);
            
            gexp_n = size(Uhat_clust,1);
            for c_idx = 1:cluster_n
                obj.cluster_sz(c_idx) = nnz(obj.gexp_cluster == c_idx);
                obj.cluster_w(c_idx) = sum((obj.gexp_cluster == c_idx).*[.5;.5*repmat(1/(gexp_n-1),gexp_n-1,1)]);
            end

            % create debug plots
            if obj.params.verbose
                obj.plot_cluster(Uhat_clust)
            end
        end
        
        
        function [obj, samp_data_out] = get_lexp_samp(obj, samp_n) 
            % Generate cluster_n samples with size samp_n of (hour, gexp)
            % coordinates and the w vector that will be used for the TSTNEP
            % optimization 
           

            %initialize storage            
            clust_lin_samp = zeros(obj.cluster_n, samp_n);
            clust_hr = zeros(obj.cluster_n, samp_n);
            clust_gexp = zeros(obj.cluster_n, samp_n);
            w = zeros(obj.cluster_n, samp_n);
            w_pos = zeros(obj.cluster_n, samp_n);
            clust_samp_mean = zeros(obj.cluster_n, samp_n);
            
            clust_mean = zeros(obj.cluster_n,1);
            
            % Cluster id will be useful when clusters have different samp_n
            clust_id = [];
            
            % calculations that can be performed out of the loop
            lf_n = size(obj.Uhat_gexp,2);
            Zgn = obj.unfold(2, obj.Z);
            Zgn_mean = mean(Zgn);
            G_apx = obj.Ugn(:,1:lf_n)'*(Zgn-Zgn_mean);
            
            % for each cluster take a D-opt sample of hours and gen exp plans
            for c_idx = 1:obj.cluster_n
                % get gexp in this cluster
                clust = (obj.gexp_cluster == c_idx);
                
                % make approximate A matrix for cluster of gen expansions
                Z_apx = obj.Uhat_gexp(clust,:)*G_apx + Zgn_mean;
                fold_sz = obj.size_z;
                fold_sz(2) = sum(clust);
                Z_apx = (obj.unfold(3,obj.fold(2, fold_sz, Z_apx)));
                Ztl_mean = mean(Z_apx);
                A_apx = Z_apx - Ztl_mean;

                clust_mean(c_idx) = mean(Ztl_mean);
                
                % get candidate design points for cluster
                %Xc_nosvd = (obj.Utl(:,1:(samp_n-1))'*Xc_cand)';
                [~,S,V] = svd(A_apx,'econ');
                cand_lf = min((samp_n-1), size(S,1)-1);
                Xc_cand = (S(1:cand_lf,1:cand_lf)*V(:,1:cand_lf)')';
                
                
                % find best sample for cluster
                [Xc, clust_lin_samp(c_idx,:)] = obj.D_opt_samp(Xc_cand, samp_n, 5);
                
                % record samp_n for cluster and samp mean for cluster
                clust_id = [clust_id; c_idx*ones(samp_n,1)];
                clust_samp_mean(c_idx,:) = Ztl_mean(clust_lin_samp(c_idx,:));
                
                % find hr and gexp for each sampled unit
                clust_hr(c_idx,:) = mod(clust_lin_samp(c_idx,:)-1,obj.size_z(1))+1;
                gexp_list = find(clust);
                clust_gexp(c_idx,:) = gexp_list(ceil(clust_lin_samp(c_idx,:)./obj.size_z(1)));
                
                % find w with hat matrix
                Hc = Xc_cand*((Xc'*Xc)\(Xc'));
                w(c_idx, :) = sum(Hc)./(obj.size_z(1)*obj.cluster_sz(c_idx));
                                
                % find w pos with gradient decent
                exp_cost = sum(A_apx,2)./(obj.size_z(1)*obj.cluster_sz(c_idx));
                A_apx_samp = A_apx(:,clust_lin_samp(c_idx,:));
                w_pos(c_idx,:) = obj.grad_dec_wpos(exp_cost, A_apx_samp, w(c_idx, :));
 
                
                if obj.params.verbose
                %fun = @(x)sum((A_apx_samp*x-exp_cost).^2);
                %x = fmincon(fun,abs(clust_w(c_idx, :)*10)',[],[],[],[],abs(clust_w(c_idx, :)*.01)',[]);
                    hold off
                    plot(exp_cost + clust_mean(c_idx))
                    hold on
                    plot(A_apx_samp*w_pos(c_idx, :)'+ clust_mean(c_idx))
                    plot(A_apx_samp*w(c_idx, :)'+ clust_mean(c_idx))
                    %plot(A_apx_samp*x)
                end
            end         
            
            % output data
            obj.lexp_samp_n = samp_n;
            obj.lexp_samp_id = clust_lin_samp;
            obj.lexp_hr_samp = clust_hr;
            obj.lexp_gexp_samp = clust_gexp;
            obj.lexp_clust_samp = clust_id;
            obj.lexp_samp_mean = clust_samp_mean;
            obj.cluster_mean = clust_mean;
            obj.lexp_w = w;
            obj.lexp_wpos = w_pos;
            
            samp_data_out.lexp_samp_n = samp_n;
            samp_data_out.lexp_samp_id = clust_lin_samp;
            samp_data_out.lexp_hr_samp = clust_hr;
            samp_data_out.lexp_gexp_samp = clust_gexp;
            samp_data_out.lexp_clust_samp = clust_id;
            samp_data_out.lexp_samp_mean = clust_samp_mean;
            samp_data_out.cluster_mean = clust_mean;
            samp_data_out.lexp_w = w;
            samp_data_out.lexp_wpos = w_pos;
            samp_data_out.cluster_w = obj.cluster_w;
        end
        
        
        function params_out = make_tstnep_params(obj, all_gen, sn)
            % Creates the parameter file for the TSTNEP optimization. If a
            % scenario 'sn' is specified, the parameters will be for the
            % TSTNEP model that assumes that only that cluster will occur
            % in the second stage
            
            if nargin < 3          
                params_out.hr_run_map = obj.lexp_hr_samp';
                params_out.hr_run_map = params_out.hr_run_map(:);

                params_out.gexp_run_map = obj.lexp_gexp_samp';
                params_out.gexp_run_map = params_out.gexp_run_map(:);

                params_out.clust_run_map = obj.lexp_clust_samp;

                params_out.cluster_n = obj.cluster_n;
                params_out.baseline_clust = obj.gexp_cluster(1);

                
                gen_list = obj.lexp_gexp_samp';
                params_out.gexp_gens = all_gen(gen_list(:),:);


                params_out.cluster_w = obj.cluster_w;

                params_out.run_w_pos = obj.lexp_wpos';
                params_out.run_w_pos = params_out.run_w_pos(:);

                params_out.run_w = obj.lexp_w';
                params_out.run_w = params_out.run_w(:);


                params_out.run_mean_center = obj.lexp_samp_mean';
                params_out.run_mean_center = params_out.run_mean_center(:);            

                params_out.cluster_mean = obj.cluster_mean;

                params_out.gexp_cluster_save = obj.gexp_cluster;
            else
                params_out.hr_run_map = obj.lexp_hr_samp([1,sn],:)';
                params_out.hr_run_map = params_out.hr_run_map(:);

                params_out.gexp_run_map = obj.lexp_gexp_samp([1,sn],:)';
                params_out.gexp_run_map = params_out.gexp_run_map(:);

                params_out.clust_run_map = repelem([1;2], obj.lexp_samp_n, 1);

                params_out.cluster_n = 2;
                params_out.baseline_clust = obj.gexp_cluster(1);


                gen_list = obj.lexp_gexp_samp([1,sn],:)';
                params_out.gexp_gens = all_gen(gen_list(:),:);


                params_out.cluster_w = [.5;.5];

                params_out.run_w_pos = obj.lexp_wpos([1,sn],:)';
                params_out.run_w_pos = params_out.run_w_pos(:);

                params_out.run_w = obj.lexp_w([1,sn],:)';
                params_out.run_w = params_out.run_w(:);

                params_out.run_mean_center = obj.lexp_samp_mean([1,sn],:)';
                params_out.run_mean_center = params_out.run_mean_center(:);            

                params_out.cluster_mean = obj.cluster_mean([1,sn]);

                params_out.gexp_cluster_save = obj.gexp_cluster;
                
            end
        end
        
        
        function plot_Uhat_gexp(obj, beta_maker)
            % Creates test plots to visually check Uhat extimation
            
            Zgn = obj.unfold(2, obj.Z);
            Zgn_mean = mean(Zgn);
            
            test_samp = Zgn- Zgn_mean;
            test_samp = test_samp(:,obj.gexp_samp_id);
            test_beta = (beta_maker*test_samp')'; 
            
            
            scatter(obj.Ugn(:,1), obj.Ugn(:,2));
            hold on 
            scatter(obj.Uhat_gexp(:,1),obj.Uhat_gexp(:,2)); 
            scatter(test_beta(:,1), test_beta(:,2));
        end
        
        
        function plot_cluster(obj, Uhat_clust)
            % Creates test plots to visually check clustering performance
            
            clust = find(obj.gexp_cluster == 1);
            scatter3(Uhat_clust(clust,1), Uhat_clust(clust,2), Uhat_clust(clust,3));
            hold on
            xlim([-300000000 300000000])
            ylim([-300000000 300000000])
            zlim([-300000000 300000000])
            
            for c_idx = 2:obj.cluster_n
                clust = find(obj.gexp_cluster == c_idx);
                scatter3(Uhat_clust(clust,1), Uhat_clust(clust,2), Uhat_clust(clust,3));
            end
            grid on  
        end
        

    end
    
    
    methods
        % These support methods are candidates to move out of the LUCA
        % object into more general functions
        
        function [unfolded, sz]  = unfold(obj, dim, data_in)
            % Unfolds a high order array into a 2D matrix with the 'dim'
            % dimension down the rows and all other dimensions folded into
            % the columns. Outputs the unfolded matrix and the size 'sz' of
            % the original array for refolding.
            
            if nargin  < 3
                x = obj.Z;
            else
                x = data_in;
            end
            
            % get order of data array
            sz = size(x);
            ord = length(sz);
            if dim > ord
                fprintf('folded dimention does not exist')
            end
            
            % orient data to folded dimension
            dim_order = circshift(1:ord, -(dim-1));
            x = permute(x,dim_order);

            % unfold
            r = size(x, 1);
            n = numel(x);
            unfolded = reshape(x, r, n/r);  

        end
        
        function folded = fold(~, dim, sz, data_in)
            % Takes a 2D matrix and folds it into a higher order array with
            % dimensions given by 'sz' and the rows moved to the 'dim'
            % demension
            
            % get the proper ordering for the dimensions
            % **NOTE** sz is the size of the dimensions for the destination
            % array so we need to order it for the array that will come out
            % of the reshape operation then dim_order it back to the new
            % array
            
            ord = length(sz);
            sz_order = circshift(1:ord, - (dim-1));
            dim_order = circshift(1:ord, (dim-1));
            ordered_sz = sz(sz_order);
            
            if dim > ord
                fprintf('folded dimention does not exist')
            elseif prod(sz) ~= numel(data_in)
                fprintf('Array does not fit into the given demensions');
            end
            
            x = data_in;
            x = reshape(x, ordered_sz);

            folded = permute(x, dim_order);
                 
        end
        
        function [X, X_id] = D_opt_samp(~, X_cand, samp_n, tries)
            % Takes candidate design points X_cand and number of samples
            % needed and returns a near D optimal design matrix X and ids
            % that correspond to the candidate matrix 
            
            if nargin < 4
                tries = 20;
            end
            
            % take initial sample
            samp = candexch(X_cand, samp_n, 'tries', tries,'display', 'on');

            % if needed find additional samples
            samp = unique(samp);
            while length(samp) < samp_n
                [X_new, X_id_map] = setdiff(X_cand, X_cand(samp,:), 'rows');
                new_samps = candexch(X_new, samp_n - length(samp), 'start', X_cand(samp,:), 'tries', tries, 'display', 'off');
                samp = unique([samp; X_id_map(new_samps)]);
            end
            
            X_id = samp;
            X = X_cand(samp,:);
%todo : evaluate D and G optimality of sample  

        end
        
        function [w_pos, fail_indicator] = grad_dec_wpos(~, target, X, w_start)
            % finds positive weights w_pos that best approximate the 
            % 'target' vector with input values X s.t. target = X*w_pos and
            % w_pos > 0
                     

            % calculate positive weights with penalized gradient descent
            reg_w = w_start;
            neg_w = find(w_start < 0.001);
            pos_w = find(w_start > 0.001);
            [row_n,~] = size(X);
            count = 0;

            % continue until no negative weights exist or convergence error 
            while (~isempty(neg_w)) && (count < 750)
                %penalize negative weights
                reg_w(neg_w) = reg_w(neg_w) + abs(w_start(neg_w)*.02*(count/100));

                % adjust positive weights to accommodate inflating negative values
                for d_idx = 1:100
                    err = X*reg_w' - target;
                    row_order = randperm(row_n);
                    alpha = max(abs(reg_w))/(max(count,100)*max(err)*max(sum(abs(X(:,pos_w)),2)));
                    for r_idx = row_order
                        reg_w(pos_w) = reg_w(pos_w) - alpha*err(r_idx)*X(r_idx,pos_w);
                    end
                end

                % check for positive and negative weights
                neg_w = find(reg_w < 0.001);
                pos_w = find(reg_w > 0.001);
                count = count+1;
            end

            % assign resulting positive weights
            if sum(sign(reg_w))~= length(reg_w)
                fail_indicator = 1;
                w_pos = max(0,w_start);
            else
                fail_indicator = 0;
                w_pos = reg_w;
            end
        end
    end
end