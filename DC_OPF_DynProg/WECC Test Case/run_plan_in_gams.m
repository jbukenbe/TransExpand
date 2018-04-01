function problem = run_plan_in_gams(problem)
% this function takes the problem input from matlab and formats it to run
% the opf problem in gams and output the cost from each scenario

%History            
%Version    Date        Who     Summary
%1          03/09/2018  JesseB  Initial Version


%% Initialize Data
candidate_plan = problem.lines;
%scen_n = problem.params.scen.n;
scen_n = problem.scen;
bus_n = problem.params.bus.n;
gen_n = problem.params.gen.n;
non_dispatch_gen = problem.non_dispatch;
scen_offset = 1000;
scen_r = scen_offset:(scen_offset+scen_n-1);

%% Format Sets for GAMS
% set string creaton function
get_uel = (@(s,v) strcat(s,strsplit(num2str(v))));

% subset of scenarios to run
subs.name='subs';
subs.uels = get_uel('s',scen_r);

% subset of lines to include in candidate plan
lines.name = 'lines';
lines.uels = get_uel('l',1:candidate_plan);

% parameter of loads at each bus and scenario pair
pload.name = 'pload';
pload.type = 'parameter';
pload.form = 'full';
pload.uels={get_uel('i', 1:bus_n), subs.uels};
load_ag = matfile('area_load_data.mat');
bus_disag = matfile('bus_area_data.mat');
pload.val = bus_disag.bus_area*load_ag.area_load(scen_r,:)';

% parameter of renewable generation at non-dispatchable scenario pair
pRenew.name = 'pRenew';
pRenew.type = 'parameter';
pRenew.form = 'full';
pRenew.uels = {get_uel('g',non_dispatch_gen), subs.uels};
renew_gen_profile = matfile('renew_dispatch_data.mat');
renew_gen_idx = load('renew_gen_profile_index.mat');
r_gen_idx = renew_gen_idx.renew_idx(renew_gen_idx.renew_idx ~=0);
pRenew.val = renew_gen_profile.renew_data(:,scen_r);
pRenew.val = max(0,pRenew.val(r_gen_idx,:));
g_diff = length(non_dispatch_gen)-length(r_gen_idx);
pRenew.val = [zeros(g_diff, scen_n);pRenew.val];

% Parameter for gen costs over scenarios
pVarCost.name = 'pVarCost';
pVarCost.type = 'parameter';
pVarCost.form = 'full';
pVarCost.uels = {get_uel('g',1:gen_n), subs.uels};
load('gen_cost_data.mat');
monthly_gen_fuel_cost = gen_cost.heatrate.*gen_cost.fuel_cost;
scen_month = ceil((scen_r)./720);
pVarCost.val = monthly_gen_fuel_cost(:,scen_month);
pVarCost.val = pVarCost.val + repmat(gen_cost.vom,1,scen_n);

%% Run GAMS Optimiazation
filename = ['MtoG', '.gdx'];
wgdx(filename, subs, lines, pload, pRenew, pVarCost);
system ('gams "OPF_file" lo=3');

%% Get GAMS Output
gams_out.name='scen_cost';
gams_out.form='full';
output = rgdx('results',gams_out);

problem.scen_op_cost = output.val;
    
end


