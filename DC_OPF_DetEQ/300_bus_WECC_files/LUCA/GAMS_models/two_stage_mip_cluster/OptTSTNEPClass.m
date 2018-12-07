classdef OptTSTNEPClass
% This class is used to handle the Two Stage Transmission Network Expansion
% Planning (TSTNEP) problem. An instance of this class will store the
% needed data and execute optimizations when needed. 
        
%History            
%Version    Date        Who     Summary
%1          11/13/2018 JesseB   Created (not working) 
%2          11/29/2018 JesseB   Working version

    properties
        tstnep_params        
        opf_params
        opt_params
        
        obj_val
        opt_decision
        run_time   
    end
    
    methods
        function obj = OptTSTNEPClass(tstnep_params, opf_params, opt_params)
            % Instance Creation
            if nargin < 3
                opf_params = {};            
                if nargin < 4
                    opt_params = {};
                end
            end
            obj.opf_params = opf_params;
            obj.opt_params = opt_params;
            obj.tstnep_params = tstnep_params;         
        end
        
        
        function obj = set.opf_params(obj, opf_params)
            % Set model parameters to defaults if none are specified for
            % each parameter
            
            opf_default = struct(   'load_growth', 1.2,...
                                    'hydro_scale', .25);
                                
                                
            if isempty(opf_params)
                obj.opf_params = opf_default;
            else
                fields = fieldnames(opf_default);
                for i = 1:length(fields)
                    if ~isfield(opf_params,fields{i})
                        opf_params.(fields{i}) = opf_default.(fields{i});
                    end
                end
                obj.opf_params = opf_params;
            end
        end
        
        
        function obj = set.opt_params(obj, opt_params)
            % Set CPLEX parameters to defaults if none are specified for
            % each parameter
            
            opt_default = struct(   'lpmethod', 6,...
                                    'SolnPoolPop', 1,...
                                    'mipemphasis', 0,...
                                    'mipstart', 1,...
                                    'startalg', 4,...
                                    'subalg', 0,...
                                    'epper', .0000001);
                                
                                
            if isempty(opt_params)
                obj.opt_params = opt_default;
            else
                fields = fieldnames(opt_default);
                for i = 1:length(fields)
                    if ~isfield(opt_params,fields{i})
                        opt_params.(fields{i}) = opf_default.(fields{i});
                    end
                end
                obj.opt_params = opt_params;
            end
        end  
        
        
        function [obj, opt_out] = solve(obj)
            % Run the Two Stage Transmission Network Expansion Planning
            % Problem with clusters for the given parameters 
            
            %% Load Data
            gen_data = load('gen_data.mat');  
            load bus_area_data.mat bus_area;

            %% Initialize Data
            % **NOTE** currently 'scen' refers to generation expansion scenarios
            % and sub refers to the subset of scen that is run in the model
            
            [hour_list,~,hour_idx] = unique(obj.tstnep_params.hr_run_map);
            [sub_list,~,subs_idx] = unique(obj.tstnep_params.gexp_run_map);
                        
            bus_n = size(bus_area,1);
            gen_n = size(gen_data.built,1);
            cluster_n = obj.tstnep_params.cluster_n;
            hr_n = length(hour_list);
            scen_n = length(sub_list);
            
            load_growth = obj.opf_params.load_growth*ones(1,scen_n);
            
            %% Format Sets for GAMS
            % set string creaton function
            get_uel = (@(s,v) strcat(s,strsplit(num2str(v))));

            % GAMS set definitions
            gen.uels = get_uel('g',1:gen_n);
            bus.uels = get_uel('i', 1:bus_n);


            % scenarios in model
            s.name= 's';
            s.uels = get_uel('s',1:scen_n);


            % subset of scenarios to run
            subs.name = 'subs';
            subs.uels = get_uel('s',sub_list');


            % subset of hours in model
            subh.name = 'subh';
            subh.uels = get_uel('h',hour_list');


            % subset of w clusters
            c.name = 'c';
            c.uels = get_uel('c',1:cluster_n);


            % network of w clusters
            cc.name = 'cc';
            cc.uels = {c.uels,c.uels};
            base_c = obj.tstnep_params.baseline_clust;
            cc.val = [base_c*ones(cluster_n-1,1), setdiff((1:cluster_n)',base_c)];


            % mapping of scenarios to w clusters
            csh.name = 'csh';
            csh.uels = {c.uels, subs.uels, subh.uels};
            csh.val = [obj.tstnep_params.clust_run_map, subs_idx, hour_idx];


            % mapping of generators included in each scenario
            sg.name = 'sg';
            sg.uels = {subs.uels, gen.uels};
            load gen_map_data.mat gen_cand gen_exist; 
            sg.val = [];
            for sub_idx = 1:length(sub_list)
                new_gens = obj.tstnep_params.gexp_gens(sub_idx,:)'.*gen_cand;
                new_gens(new_gens == 0) = [];
                scen_gens = sort([gen_exist;new_gens]);
                this_scen = sub_idx*ones(length(scen_gens),1);
                sg.val = [sg.val; this_scen, scen_gens];
            end
            clear gen_cand gen_exist


            % parameter of cluster weights
            pc_w.name = 'pc_w';
            pc_w.type = 'parameter';
            pc_w.form = 'full';
            pc_w.uels = c.uels;
            pc_w.val = obj.tstnep_params.cluster_w;


            % parameter of hour weights for subh
            psub_w.name = 'psub_w';
            psub_w.type = 'parameter';
            psub_w.form = 'sparse';
            psub_w.uels = csh.uels;
            psub_w.val = [csh.val, obj.tstnep_params.run_w_pos];


            % parameter of hourly means
            ph_means.name = 'ph_means';
            ph_means.type = 'parameter';
            ph_means.form = 'sparse';
            ph_means.uels = csh.uels;
            ph_means.val = [csh.val, obj.tstnep_params.run_mean_center];


            % parameter of cluster means
            pcluster_means.name = 'pcluster_means';
            pcluster_means.type = 'parameter';
            pcluster_means.form = 'full';
            pcluster_means.uels = c.uels;
            pcluster_means.val = obj.tstnep_params.cluster_mean;


            % parameter of loads at each bus, time stage, scenario, and hour
            pload.name = 'pload';
            pload.type = 'parameter';
            pload.form = 'full';
            load area_load_data.mat area_load;
            pload.uels={bus.uels, subs.uels, subh.uels};
            pload.val = bus_area*area_load(hour_list,:)';
            pload.val = bsxfun(@times,repmat(pload.val,1,1,scen_n),permute(load_growth,[3 1 2]));
            pload.val = permute(pload.val,[1,3,2]);
            clear area_load bus_area


            % parameter of renewable generation at non-dispatchable scenario pair
            % combine renewable profile with generic hydro id
            renew_gen_idx = gen_data.renew_id + gen_data.generic_hydro;
            non_dispatch_gen = find(renew_gen_idx);
            renew_gen_idx = renew_gen_idx(non_dispatch_gen);

            % subset of non-disptachable generation
            nd.name = 'nd';
            nd.uels = get_uel('g',non_dispatch_gen');

            pRenew.name = 'pRenew';
            pRenew.type = 'parameter';
            pRenew.form = 'full';
            pRenew.uels = {nd.uels, subs.uels, subh.uels};  
            
            load renew_dispatch_data.mat renew_data;
            renew_data(81,:) = renew_data(81,:).*obj.opf_params.hydro_scale;
            pRenew.val = renew_data(:,hour_list);
            clear renew_data
            pRenew.val = max(0,pRenew.val(renew_gen_idx,:));
            pRenew.val = permute(repmat(pRenew.val,1,1,scen_n),[1, 3, 2]);


            % Parameter for gen costs over scenarios
            pVarCost.name = 'pVarCost';
            pVarCost.type = 'parameter';
            pVarCost.form = 'full';
            pVarCost.uels = {gen.uels, subs.uels, subh.uels};
            hour_month = ceil((hour_list)./744);
            fuel_cost_temp = gen_data.fuel_cost;
            gen_fuel_cost = fuel_cost_temp(gen_data.fuel_id, hour_month);
            clear fuel_cost_temp
            pVarCost.val = gen_data.heatrate.*gen_fuel_cost;
            pVarCost.val = pVarCost.val + repmat(gen_data.vom,1,hr_n); 
            pVarCost.val = permute(repmat(pVarCost.val,1,1,scen_n),[1, 3, 2]);


            % generator min
            pgen_min.name = 'pgen_min';
            pgen_min.type = 'parameter';
            pgen_min.form = 'full';
            pgen_min.uels = gen.uels;
            pgen_min.val = gen_data.pmin;


            % generator max
            pgen_max.name = 'pgen_max';
            pgen_max.type = 'parameter';
            pgen_max.form = 'full';
            pgen_max.uels = gen.uels;
            pgen_max.val = gen_data.pmax;


            %% Run GAMS Optimiazation
            filename = ['MtoG', '.gdx'];
            wgdx(filename, c, subs, subh, nd, sg, cc, csh, pc_w, psub_w, ph_means, pcluster_means , pload, pRenew, pVarCost, pgen_min, pgen_max);
            tic;
            system ('gams "clust_tep2_MIP" lo=3');
            obj.run_time = toc;

            %% Get GAMS Output
            % scenario costs in optimal plan
            scen_cost.name='scen_cost';
            scen_cost.form='sparse';
            scen_cost.uels = csh.uels;
            scen_output = rgdx('results',scen_cost);

            % optimal plan total cost
            tot_cost.name='tot_cost';
            tot_cost.form='full';
            tot_output = rgdx('results',tot_cost);

            % lines built in optimal plan
            lines.uels = get_uel('l',655:705);
            lines_built.name='lines_built';
            lines_built.form='sparse';
            lines_built.uels = {lines.uels, bus.uels, bus.uels, c.uels};
            line_output = rgdx('results',lines_built);

            lb.name = 'lb';
            lb.form = 'full';
            lb_output = rgdx('results',lb);
            
            opt_out.scen_cost = scen_output.val;
            opt_out.total_cost = tot_output.val;
            opt_out.line_out = line_output.val;
            opt_out.lb = lb_output.val;
            
            obj.obj_val = tot_output.val;
            obj.opt_decision = line_output.val;
            obj.run_time = toc;
        end
    end
  
end