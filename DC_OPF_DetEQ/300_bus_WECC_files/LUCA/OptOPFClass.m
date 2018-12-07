classdef OptOPFClass
% This class is used to handle the OPF optimizations required by the
% Sampling and Approximation class of algorithms. Objects of this class
% type will store the needed data and execute optimizations
        
%History            
%Version    Date        Who     Summary
%1          11/13/2018  JesseB  Created (not working)
%2          11/16/2018  JesseB  Initial Working version
%
%
%   TODO:   Better writing of cplex parameters
%           Make 'run_all' methods efficient and general to disjoint hrs
%           

%g = gams(struct('gams','C:/GAMS/win64/25.1/gams.exe'))
    
    properties
        opf_params
        opt_params
        
        lexp_list
        gexp_list
        hr_list
        hr_n
        lexp_n
        gexp_n
        
        crossruns
        run_n
        
        cand_lexp
        cand_gexp
        
        opt_val
        opt_decision
        run_time
    end
    
    methods
        function obj = OptOPFClass(hr_in, line_in, gen_in, opf_params, opt_params)
            % Instance creation method
            % note that if parameters are not set, defaults are used
            if nargin < 4
                opf_params = {};            
                if nargin < 5
                    opt_params = {};
                end
            end
            obj.opf_params = opf_params;
            obj.opt_params = opt_params;
            
            
            % Set run length data
            obj.hr_n = length(hr_in);
            obj.lexp_n = size(line_in, 1);
            obj.gexp_n = size(gen_in, 1);
            
            if (obj.lexp_n ~= obj.hr_n) || (obj.gexp_n ~= obj.hr_n)
                obj.crossruns = 1;
                obj.run_n = obj.hr_n*obj.lexp_n*obj.gexp_n;
            else
                obj.crossruns = 0;
                obj.run_n = length(hr_in);
            end
            obj.lexp_list = line_in;
            obj.gexp_list = gen_in;
            obj.hr_list = hr_in;    
        end
        
        
        function obj = set.opf_params(obj, opf_params)
            % Set model parameters to defaults if none are specified for
            % each parameter
            
            opf_default = struct(   'load_growth', 1.2,...
                                    'hydro_scale', .25,...
                                    'group_size', 12);
                                
                                
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
        
        
        function cand_lexp = get_cand_lexp(obj, lexp_id)
            new_lines = obj.lexp_list(lexp_id,:).*(655:705);
            new_lines(new_lines == 0) = [];
            cand_lexp = [1:654, new_lines];    
        end
        
        
        function cand_gexp = get_cand_gexp(obj, gexp_id)
            load gen_map_data.mat gen_cand gen_exist; 
            new_gens = obj.gexp_list(gexp_id,:).*gen_cand';
            new_gens(new_gens == 0) = [];
            cand_gexp = sort([gen_exist' , new_gens]);
        end
        
        
        function [obj, opt_val] = solve_plan(obj, lexp_id, gexp_id)
            % solves the opf for every hour of the plan specified by the
            % line expansion id 'lexp_id' and generation expansion id
            % 'gexp_id'. These ID are from the order given to the object.
            
% todo make this actually read the parameters             
            opt_param = [6 2 0 1 4 0 .0000001];
            fileID = fopen('cplex.opt','w');
            formatSpec = 'lpmethod %d\nSolnPoolPop %d\nmipemphasis %d\nmipstart %d\nstartalg %d\nsubalg %d\nepper %7.7f\n';
            fprintf(fileID,formatSpec,opt_param);
            fclose(fileID);

            
            % make groups of hours
            group_size = obj.opf_params.group_size;
            plan_run_n = ceil(obj.hr_n/group_size);
            
            % loop through groups of hours
            opf_out = zeros(obj.hr_n,1);
            for r_idx = 1:plan_run_n
                if r_idx ~= plan_run_n
                    save_idx = (1+group_size*(r_idx-1)):(group_size*(r_idx));
                else
                    save_idx = (1+group_size*(r_idx-1)):obj.hr_n;                   
                end
                hr_id = obj.hr_list(save_idx);
                
                % solve opf group
                [~,opf_out(save_idx)] = obj.solve(hr_id, lexp_id, gexp_id);
            end
            % output results
            opt_val = opf_out;
        end

        
        function [obj, opt_val] = solve(obj, hr_id, lexp_id, gexp_id)
            % Solves the OPF problem for the given hour and expansion plan
            % note that the currently only one line and generation
            % expansion id is supported but multiple hr_ids can be given
            % but are classified as scenarios from legacy code.
            
            %% Load Data
            load_ag = matfile('area_load_data.mat');
            bus_area_data = matfile('bus_area_data.mat');
            renew_gen_profile = matfile('renew_dispatch_data.mat');
            gen_data = load('gen_data.mat');  
            
            %% Initialize Data
            bus_n = size(bus_area_data.bus_area,1);
            gen_n = size(gen_data.built,1);
            
            load_growth = obj.opf_params.load_growth;
            hydro_scale = obj.opf_params.hydro_scale;
                        
            gen_plan = obj.get_cand_gexp(gexp_id);            
            line_plan = obj.get_cand_lexp(lexp_id);

            scen_list = hr_id;   
            scen_n = length(scen_list);
            scen_w = ones(1,scen_n).*(1/8760);

        %% Format Sets for GAMS
            % set string creaton function
            get_uel = (@(s,v) strcat(s,strsplit(num2str(v))));

            % gen and bus sets
            gen.uels = get_uel('g',1:gen_n);
            bus.uels = get_uel('i', 1:bus_n);

            % subset of scenarios to run
            subs.name= 'subs';
            subs.uels = get_uel('s',scen_list);

            % subset of lines to include in candidate plan
            lines.name = 'lines';
            lines.uels = get_uel('l',line_plan);

            % subset of generator to run
            gens.name = 'gens';
            gens.uels = get_uel('g',gen_plan);


            %% Format Parameters for GAMS
            % parameter of scenario weights for subs
            psub_w.name = 'psub_w';
            psub_w.type = 'parameter';
            psub_w.form = 'full';
            psub_w.uels = subs.uels;
            psub_w.val = scen_w;

            % parameter of loads at each bus and scenario pair
            pload.name = 'pload';
            pload.type = 'parameter';
            pload.form = 'full';
            pload.uels={bus.uels, subs.uels};
            temp_load_ag = load_ag.area_load;
            pload.val = load_growth.*bus_area_data.bus_area*temp_load_ag(scen_list,:)';
            clear temp_load_ag

            % get variable op costs
            pVarCost.name = 'pVarCost';
            pVarCost.type = 'parameter';
            pVarCost.form = 'full';
            pVarCost.uels = {gen.uels, subs.uels};
            scen_month = ceil((scen_list)./744);
            fuel_cost_temp = gen_data.fuel_cost;
            gen_fuel_cost = fuel_cost_temp(gen_data.fuel_id, scen_month);
            clear fuel_cost_temp
            pVarCost.val = gen_data.heatrate.*gen_fuel_cost;
            pVarCost.val = pVarCost.val + repmat(gen_data.vom,1,scen_n);    


            % get renewable profiles
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
            pRenew.uels = {nd.uels, subs.uels};    
            temp_renew_data = renew_gen_profile.renew_data;
            % duct tape correction for hydro data being too large
            %temp_renew_data(81,:) = temp_renew_data(81,:).*.5;
            % better solution to hydro:
            temp_renew_data(81,:) = temp_renew_data(81,:).*hydro_scale;
            pRenew.val = temp_renew_data(:,scen_list);
            clear temp_renew_data
            pRenew.val = max(0,pRenew.val(renew_gen_idx,:));


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


            %% Write GAMS input data file and run optimization
            filename = ['MtoG', '.gdx'];
            wgdx(filename, subs, lines, gens, nd, pgen_min, pgen_max, psub_w, pload, pRenew, pVarCost);
            system ('gams "OPF_file" lo=2');


            %% Get GAMS Output
            % scenario costs 
            scen_cost.name='scen_cost';
            scen_cost.form='full';
            scen_cost.uels = {subs.uels};
            scen_output = rgdx('results',scen_cost);

            % shed generation
            p_shed_gen.name='p_shed_gen';
            p_shed_gen.form='full';
            p_shed_gen.uels = {gen.uels, subs.uels};
            shed_gen = rgdx('results', p_shed_gen);

            % power not served
            p_pns.name='p_pns';
            p_pns.form='full';
            p_pns.uels = {bus.uels, subs.uels};
            pns_output = rgdx('results',p_pns);    


            opt_val = scen_output.val;
            
        end
    end
end