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

%% Format Sets for GAMS
% set string creaton function
get_uel = (@(s,v) strcat(s,strsplit(num2str(v))));

% subset of scenarios to run
subs.name='subs';
subs.uels = get_uel('s',1:scen_n);

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
pload.val = bus_disag.bus_area*load_ag.area_load(1:scen_n,:)';

% parameter of renewable generation at non-dispatchable scenario pair
pRenew.name = 'pRenew';
pRenew.type = 'parameter';
pRenew.form = 'full';
%pRenew.uels = {get_uel('g',non_dispatch_gen), subs.uels};
%renew_gen = matfile('renew_gen_data.mat');
%pRenew.val = renew_gen.dispatch(:,1:scen_n);


% Parameter for gen costs over scenarios
pVarCost.name = 'pVarCost';
pVarCost.type = 'parameter';
pVarCost.form = 'full';
pVarCost.uels = {get_uel('g',1:gen_n), subs.uels};
load('gen_cost_data.mat');
monthly_gen_fuel_cost = gen_cost.heatrate.*gen_cost.fuel_cost;
scen_month = ceil((1:scen_n)./30);
pVarCost.val = monthly_gen_fuel_cost(:,scen_month);
pVarCost.val = pVarCost.val + repmat(gen_cost.vom,1,scen_n);

%% Run GAMS Optimiazation
filename = ['MtoG', '.gdx'];
wgdx(filename, subs, lines, pload);
system ('gams "OPF_file" lo=2');

%% Get GAMS Output
gams_out.name='scen_cost';
gams_out.form='full';
output = rgdx('results',gams_out);

problem.scen_op_cost = output.val;
    
end


