function problem = run_tep_test_gams(problem)
% this function takes the problem input for the transmission expansion
% problem and formulates it into a GAMS MIP that solves the single stage
% stochastic TEP. A scenario subset and their weights are needed for the
% optimization to be compatable with compressed models.


%History            
%Version    Date        Who     Summary
%1          06/05/2018  JesseB  Adapted from run_opf_gams
%2          06/11/2018  JesseB  Added mean parameters for approximation
%2          07/09/2018  JesseB  Bender and MIP compatable

%g = gams(struct('gams','C:/GAMS/win64/25.1/gams.exe'))

%% Initialize Data
bus_n = problem.bus_n;
gen_n = problem.gen_n;
problem.non_dispatch_gen = 332:980;
non_dispatch_gen = problem.non_dispatch_gen;
scen_w = problem.scen_w;
scen_list = problem.scen_list;
load_growth = problem.load_growth;
scen_n = length(scen_list);

%% Format Sets for GAMS
% set string creaton function
get_uel = (@(s,v) strcat(s,strsplit(num2str(v))));

% gen and bus sets
gen.uels = get_uel('g',1:gen_n);
bus.uels = get_uel('i', 1:bus_n);

% subset of scenarios to run
subs.name= 'subs';
subs.uels = get_uel('s',scen_list);


%% Format Parameters for GAMS
% parameter of scenario weights for subs
psub_w.name = 'psub_w';
psub_w.type = 'parameter';
psub_w.form = 'full';
psub_w.uels = subs.uels;
psub_w.val = scen_w;

% k means model has stochastic data saved as cluster centroids in problem structure
if problem.samp_method == "KM"
% parameter of loads at each bus and scenario pair
    pload.name = 'pload';
    pload.type = 'parameter';
    pload.form = 'full';
    pload.uels={bus.uels, subs.uels};
    pload.val = problem.scen_load;

% parameter of renewable generation at non-dispatchable scenario pair
    pRenew.name = 'pRenew';
    pRenew.type = 'parameter';
    pRenew.form = 'full';
    pRenew.uels = {get_uel('g',non_dispatch_gen), subs.uels};
    pRenew.val = problem.scen_renew;
    
% parameter for gen costs over scenarios
    pVarCost.name = 'pVarCost';
    pVarCost.type = 'parameter';
    pVarCost.form = 'full';
    pVarCost.uels = {gen.uels, subs.uels};
    pVarCost.val = problem.scen_VarCost;
    
% if k means is not used, calculate the scenario values from data
else    
% preload model data files
    load_ag = matfile('area_load_data.mat');
    bus_disag = matfile('bus_area_data.mat');
    renew_gen_profile = matfile('renew_dispatch_data.mat');
    renew_gen_idx = load('renew_gen_profile_index.mat');
    gen_data = load('gen_cost_data.mat');     
    
% parameters exclusive to latent factor approximation model
    if problem.samp_method == "LF"
        % problem global mean estimate
        global_mean.name = 'global_mean';
        global_mean.type = 'parameter';
        global_mean.form = 'full';
        global_mean.uels = {'m1'};
        global_mean.val = problem.global_mean;

        % scenario mean estimates
        psub_mean.name= 'psub_mean';
        psub_mean.type = 'parameter';
        psub_mean.form = 'full';
        psub_mean.uels = subs.uels;
        psub_mean.val = problem.psub_m;
    end

% parameter of loads at each bus and scenario pair
    pload.name = 'pload';
    pload.type = 'parameter';
    pload.form = 'full';
    pload.uels={bus.uels, subs.uels};
    temp_load_ag = load_ag.area_load;
    pload.val = load_growth.*bus_disag.bus_area*temp_load_ag(scen_list,:)';
    clear temp_load_ag

% parameter of renewable generation at non-dispatchable scenario pair
    pRenew.name = 'pRenew';
    pRenew.type = 'parameter';
    pRenew.form = 'full';
    pRenew.uels = {get_uel('g',non_dispatch_gen), subs.uels};
    r_gen_idx = renew_gen_idx.renew_idx(renew_gen_idx.renew_idx ~=0);
    temp_renew_data = renew_gen_profile.renew_data;
    pRenew.val = temp_renew_data(:,scen_list);
    clear temp_renew_data
    pRenew.val = max(0,pRenew.val(r_gen_idx,:));
    g_diff = length(non_dispatch_gen)-length(r_gen_idx);
    pRenew.val = [zeros(g_diff, scen_n);pRenew.val];
    % duct tape correction for hydro data being too large
    pRenew.val(1:107,:) = pRenew.val(1:107,:).*.5;

% parameter for gen costs over scenarios
    pVarCost.name = 'pVarCost';
    pVarCost.type = 'parameter';
    pVarCost.form = 'full';
    pVarCost.uels = {gen.uels, subs.uels};
    monthly_gen_fuel_cost = gen_data.gen_cost.heatrate.*gen_data.gen_cost.fuel_cost;
    scen_month = ceil((scen_list)./744);
    pVarCost.val = monthly_gen_fuel_cost(:,scen_month);
    pVarCost.val = pVarCost.val + repmat(gen_data.gen_cost.vom,1,scen_n);
end


%% Write GAMS input data file and run optimization
filename = ['MtoG', '.gdx'];
tic
if problem.samp_method == "LF"
    wgdx(filename, subs, psub_w, psub_mean, global_mean, pload, pRenew, pVarCost);
    if problem.opt_mod == "bender"
        system ('gams "lf_tep_bender" lo=2');
    elseif problem.opt_mod == "MIP"
        system ('gams "lf_tep_MIP" lo=2');
    end
else
    wgdx(filename, subs, psub_w, pload, pRenew, pVarCost);
    if problem.opt_mod == "bender"
        system ('gams "tep_bender" lo=2');
    elseif problem.opt_mod == "MIP"
        system ('gams "tep_MIP" lo=2');
    end
end
problem.run_time = toc;

%% Get GAMS Output
% scenario costs in optimal plan
scen_cost.name='scen_cost';
scen_cost.form='full';
scen_cost.uels = {subs.uels};
scen_output = rgdx('results',scen_cost);
problem.optimal_scen_op_cost = scen_output.val;

% lines built in optimal plan
lines.uels = get_uel('l',655:705);
lines_built.name='lines_built';
lines_built.form='sparse';
lines_built.uels = {lines.uels,bus.uels,bus.uels};
line_output = rgdx('results',lines_built);
problem.lines = line_output.val(:,1);

% optimal plan total cost
tot_cost.name='tot_cost';
tot_cost.form='full';
tot_output = rgdx('results',tot_cost);
problem.tot_cost = tot_output.val;

% optimal plan lower bound
lb.name = 'lb';
lb.form = 'full';
lb_output = rgdx('results',lb);
problem.lower_bound = lb_output.val;

%optimal plan optimality gap
problem.opt_gap = (problem.tot_cost-problem.lower_bound)/(abs(problem.lower_bound)+1);

end














