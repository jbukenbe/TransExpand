function problem = run_tep2_gams(problem)
% this function takes the problem input for the transmission expansion
% problem and formulates it into a GAMS MIP that solves the two stage
% stochastic TEP. Clustering information for the second stage decision
% nodes is needed for compressed models.


%History            
%Version    Date        Who     Summary
%1          05/23/2018  JesseB  Adapted from run_tep_gams

%TO DO:
%   update for t,s,h index
%   add data for t2 scens
%   update problem inputs for t2




%% Initialize Data
%candidate_plan = problem.candidate_plan;
%cand_n = 51;
bus_n = problem.bus_n;
gen_n = problem.gen_n;
hour_n = 3;
scen_n = 5;
time_n = 2;
non_dispatch_gen = problem.non_dispatch_gen;
%load_growth = problem.load_growth;
problem.scen.load_growth = 1.40:-(.2/(scen_n-1)):1.2;
scen_list = 1:scen_n;
hour_list = 1:hour_n;
scens_w = ones(scen_n,1)./scen_n;
scenh_w = ones(hour_n,1).*8760./hour_n;
%scen_list = problem.scen_list;
%scen_w = problem.scen_w;
scen_n = length(scen_list);

%% Format Sets for GAMS
%g = gams(struct('gams','C:/GAMS/gams.exe'))

% preload model data
load_ag = matfile('area_load_data.mat');
bus_disag = matfile('bus_area_data.mat');
renew_gen_profile = matfile('renew_dispatch_data.mat');
renew_gen_idx = load('renew_gen_profile_index.mat');
gen_data = load('gen_cost_data.mat');  

% set string creaton function
get_uel = (@(s,v) strcat(s,strsplit(num2str(v))));

% GAMS set definitions
gen.uels = get_uel('g',1:gen_n);
nd_gen.uels = get_uel('g',non_dispatch_gen);
bus.uels = get_uel('i', 1:bus_n);
time.uels = get_uel('t',1:2);

% subset of scenarios to run
subs.name= 'subs';
subs.uels = get_uel('s',scen_list);

%subset of hours in scenario
subh.name = 'subh';
subh.uels = get_uel('h',hour_list);

% parameter of scenario weights for subs
psubs_w.name = 'psubs_w';
psubs_w.type = 'parameter';
psubs_w.form = 'full';
psubs_w.uels = subs.uels;
psubs_w.val = scens_w;

% parameter of hour weights for subh
psubh_w.name = 'psubh_w';
psubh_w.type = 'parameter';
psubh_w.form = 'full';
psubh_w.uels = subh.uels;
psubh_w.val = scenh_w;

% parameter of loads at each bus, time stage, scenario, and hour
pload.name = 'pload';
pload.type = 'parameter';
pload.form = 'full';
pload.uels={bus.uels, time.uels, subs.uels, subh.uels};
temp_load_ag = load_ag.area_load;
pload.val = bus_disag.bus_area*temp_load_ag(hour_list,:)';
pload.val = bsxfun(@times,repmat(pload.val,1,1,scen_n,time_n),permute(problem.scen.load_growth,[3 4 2 1]));
pload.val = permute(pload.val,[1,4,3,2]);
clear temp_load_ag

% parameter of renewable generation at non-dispatchable scenario pair
pRenew.name = 'pRenew';
pRenew.type = 'parameter';
pRenew.form = 'full';
pRenew.uels = {nd_gen.uels, time.uels, subs.uels, subh.uels};
r_gen_idx = renew_gen_idx.renew_idx(renew_gen_idx.renew_idx ~=0);
temp_renew_data = renew_gen_profile.renew_data;
pRenew.val = temp_renew_data(:,hour_list);
clear temp_renew_data
pRenew.val = max(0,pRenew.val(r_gen_idx,:));
g_diff = length(non_dispatch_gen)-length(r_gen_idx);
pRenew.val = [zeros(g_diff, hour_n);pRenew.val];
% duct tape correction for hydro data being too large
pRenew.val(1:107,:) = pRenew.val(1:107,:).*.5;
pRenew.val = permute(repmat(pRenew.val,1,1,scen_n,time_n),[1, 4, 3, 2]);

% Parameter for gen costs over scenarios
pVarCost.name = 'pVarCost';
pVarCost.type = 'parameter';
pVarCost.form = 'full';
pVarCost.uels = {gen.uels, time.uels, subs.uels, subh.uels};
monthly_gen_fuel_cost = gen_data.gen_cost.heatrate.*gen_data.gen_cost.fuel_cost;
scen_month = ceil((hour_list)./744);
pVarCost.val = monthly_gen_fuel_cost(:,scen_month);
pVarCost.val = pVarCost.val + repmat(gen_data.gen_cost.vom,1,hour_n);
pVarCost.val = permute(repmat(pVarCost.val,1,1,scen_n,time_n),[1, 4, 3, 2]);

%% Run GAMS Optimiazation
filename = ['MtoG', '.gdx'];
wgdx(filename, subs, subh, psubs_w, psubh_w, pload, pRenew, pVarCost);
system ('gams "tep2_MIP" lo=3');


%% Get GAMS Output
% scenario costs in optimal plan
scen_cost.name='scen_cost';
scen_cost.form='full';
scen_cost.uels = {time.uels, subs.uels, subh.uels};
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
lines_built.uels = {lines.uels,bus.uels,bus.uels,time.uels, subs.uels};
line_output = rgdx('results',lines_built);
problem.lines = line_output.val(:,1);

end