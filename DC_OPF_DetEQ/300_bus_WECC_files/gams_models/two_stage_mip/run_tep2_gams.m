function problem = run_tep2_gams(problem)
% this function takes the problem input for the transmission expansion
% problem and formulates it into a GAMS MIP that solves the two stage
% stochastic TEP. Clustering information for the second stage decision
% nodes is needed for compressed models.


%History            
%Version    Date        Who     Summary
%1          05/23/2018  JesseB  Adapted from run_tep_gams


%% Initialize Data
%candidate_plan = problem.candidate_plan;
%cand_n = 51;
bus_n = problem.bus_n;
gen_n = problem.gen_n;
non_dispatch_gen = problem.non_dispatch_gen;
load_growth = problem.load_growth;
scen_n = 6;
scen_list = randperm(8736, scen_n);
scen_w = ones(scen_n,1).*8760./scen_n;
%scen_list = problem.scen_list;
%scen_w = problem.scen_w;
scen_n = length(scen_list);

%% Format Sets for GAMS
%g = gams(struct('gams','C:/GAMS/gams.exe'))

% set string creaton function
get_uel = (@(s,v) strcat(s,strsplit(num2str(v))));

% gen and bus sets
gen.uels = get_uel('g',1:gen_n);
bus.uels = get_uel('i', 1:bus_n);

% preload model data
load_ag = matfile('area_load_data.mat');
bus_disag = matfile('bus_area_data.mat');
renew_gen_profile = matfile('renew_dispatch_data.mat');
renew_gen_idx = load('renew_gen_profile_index.mat');
gen_data = load('gen_cost_data.mat');  

% subset of scenarios to run
subs.name= 'subs';
subs.uels = get_uel('s',scen_list);

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

% Parameter for gen costs over scenarios
pVarCost.name = 'pVarCost';
pVarCost.type = 'parameter';
pVarCost.form = 'full';
pVarCost.uels = {gen.uels, subs.uels};
monthly_gen_fuel_cost = gen_data.gen_cost.heatrate.*gen_data.gen_cost.fuel_cost;
scen_month = ceil((scen_list)./744);
pVarCost.val = monthly_gen_fuel_cost(:,scen_month);
pVarCost.val = pVarCost.val + repmat(gen_data.gen_cost.vom,1,scen_n);

%% Run GAMS Optimiazation
filename = ['MtoG', '.gdx'];
wgdx(filename, subs, psub_w, pload, pRenew, pVarCost);
system ('gams "transmission_exp_MIP" lo=3');


%% Get GAMS Output
% scenario costs in optimal plan
scen_cost.name='scen_cost';
scen_cost.form='full';
scen_cost.uels = {subs.uels};
scen_output = rgdx('results',scen_cost);
problem.optimal_scen_op_cost = scen_output.val;

% optimal plan total cost
tot_cost.name='tot_cost';
tot_cost.form='full';
tot_output = rgdx('results',tot_cost);
problem.tot_cost = tot_output.val;

% lines built in optimal plan
lines.uels = get_uel('l',654:705);
lines_built.name='lines_built';
lines_built.form='sparse';
lines_built.uels = {lines.uels,bus.uels,bus.uels};
line_output = rgdx('results',lines_built);
problem.lines = line_output.val(:,1);

end