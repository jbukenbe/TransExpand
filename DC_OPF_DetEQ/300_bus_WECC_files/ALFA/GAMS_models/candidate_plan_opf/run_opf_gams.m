function problem = run_opf_gams(problem)
% this function takes the problem input from matlab and formats it to run
% the opf problem in gams and output the cost from each scenario

%History            
%Version    Date        Who     Summary
%1          03/09/2018  JesseB  Initial Version
%2          04/03/2018  JesseB  Working with daily loop for speed
%3          05/18/2018  JesseB  Changed name for new file structure was run_plan_in_gams

%% Initialize Data
candidate_plan = problem.candidate_plan;
bus_n = problem.bus_n;
gen_n = problem.gen_n;
non_dispatch_gen = problem.non_dispatch_gen;
scen_n = problem.scen_n;
scen_offset = problem.scen_offset;
load_growth = problem.load_growth;

%% Format Sets for GAMS
% set string creaton function
%g = gams(struct('gams','C:/GAMS/gams.exe'))
get_uel = (@(s,v) strcat(s,strsplit(num2str(v))));

% subset of lines to include in candidate plan
lines.name = 'lines';
lines.uels = get_uel('l',candidate_plan);

gen.uels = get_uel('g',1:gen_n);
bus.uels = get_uel('i', 1:bus_n);

%% Run 6 hour loop for linnear programming efficiency
% preload model data
    load_ag = matfile('area_load_data.mat');
    bus_disag = matfile('bus_area_data.mat');
    renew_gen_profile = matfile('renew_dispatch_data.mat');
    renew_gen_idx = load('renew_gen_profile_index.mat');
    gen_data = load('gen_cost_data.mat');  
    
scen_op_cost = zeros(168,1);
shed_gen = zeros(168,1);
power_not_served = zeros(168,1);
for d_idx = 1:28
tic
    scen_r = (scen_offset+scen_n*(d_idx-1)):(scen_offset+scen_n*d_idx - 1);
    storage_r = (scen_n*(d_idx-1)+1):scen_n*d_idx;
    
% subset of scenarios to run
    subs.name='subs';
    subs.uels = get_uel('s',scen_r);    
    
% parameter of loads at each bus and scenario pair
    pload.name = 'pload';
    pload.type = 'parameter';
    pload.form = 'full';
    pload.uels={bus.uels, subs.uels};
    pload.val = load_growth.*bus_disag.bus_area*load_ag.area_load(scen_r,:)';

% parameter of renewable generation at non-dispatchable scenario pair
    pRenew.name = 'pRenew';
    pRenew.type = 'parameter';
    pRenew.form = 'full';
    pRenew.uels = {get_uel('g',non_dispatch_gen), subs.uels};
    r_gen_idx = renew_gen_idx.renew_idx(renew_gen_idx.renew_idx ~=0);
    pRenew.val = renew_gen_profile.renew_data(:,scen_r);
    pRenew.val = max(0,pRenew.val(r_gen_idx,:));
    g_diff = length(non_dispatch_gen)-length(r_gen_idx);
    pRenew.val = [zeros(g_diff, scen_n);pRenew.val];
    pRenew.val(1:107,:) = pRenew.val(1:107,:).*.5;

% Parameter for gen costs over scenarios
    pVarCost.name = 'pVarCost';
    pVarCost.type = 'parameter';
    pVarCost.form = 'full';
    pVarCost.uels = {gen.uels, subs.uels};
    monthly_gen_fuel_cost = gen_data.gen_cost.heatrate.*gen_data.gen_cost.fuel_cost;
    scen_month = ceil((scen_r)./744);
    pVarCost.val = monthly_gen_fuel_cost(:,scen_month);
    pVarCost.val = pVarCost.val + repmat(gen_data.gen_cost.vom,1,scen_n);

%% Run GAMS Optimiazation
    filename = ['MtoG', '.gdx'];
    wgdx(filename, subs, lines, pload, pRenew, pVarCost);
    system ('gams "OPF_file" lo=2');

    
%% Get GAMS Output
    scen_cost.name='scen_cost';
    scen_cost.form='full';
    scen_cost.uels = {subs.uels};
    scen_output = rgdx('results',scen_cost);
    scen_op_cost(storage_r) = scen_output.val;
    
    p_shed_gen.name='p_shed_gen';
    p_shed_gen.form='full';
    p_shed_gen.uels = {gen.uels, subs.uels};
    shed_gen_output = rgdx('results', p_shed_gen);
    shed_gen(storage_r) = sum(shed_gen_output.val);
    
    p_pns.name='p_pns';
    p_pns.form='full';
    p_pns.uels = {bus.uels, subs.uels};
    pns_output = rgdx('results',p_pns);
    power_not_served(storage_r) = sum(pns_output.val);
toc    
end
problem.scen_op_cost = scen_op_cost;
problem.shed_gen = shed_gen;
problem.pns = power_not_served;


end


