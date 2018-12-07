function problem = run_5kgen_gams(problem)
% this function takes the problem input for the WECC network and solves the
% DC OPF for the problem version with 4,992 generation units.


%History            
%Version    Date        Who     Summary
%1          07/20/2018  JesseB  Adapted from run_tep_gams

%g = gams(struct('gams','C:/GAMS/win64/25.1/gams.exe'))


%% Load Data
load_ag = matfile('area_load_data.mat');
bus_disag = matfile('bus_area_data.mat');
renew_gen_profile = matfile('renew_dispatch_data.mat');
gen_data = load('gen_data.mat');  
bus_area_data = load('bus_area_data');


%% Initialize Data
bus_n = size(bus_area_data.bus_area,1);
gen_n = size(gen_data.built,1);
load_growth = problem.load_growth;

scen_run_table = reshape(problem.scen_list,6,28);
scen_w_table = reshape(problem.scen_w,6,28);

line_plan = problem.line_plan;
gen_plan = problem.gen_plan';



%% Loop through scenarios for faster OPF solve time
problem.scen_op_cost = zeros(168,1);
for d_idx = 1:28
% Select scenarios for this run
    
    scen_w = scen_w_table(:,d_idx)./8760;
    scen_list = scen_run_table(:,d_idx)';    
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
    pload.val = load_growth.*bus_disag.bus_area*temp_load_ag(scen_list,:)';
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
    temp_renew_data(81,:) = temp_renew_data(81,:).*.5;
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
    tic

    filename = ['MtoG', '.gdx'];
    wgdx(filename, subs, lines, gens, nd, pgen_min, pgen_max, psub_w, pload, pRenew, pVarCost);
    system ('gams "OPF_file" lo=2');

    toc

    %% Get GAMS Output
    % scenario costs in optimal plan
    scen_cost.name='scen_cost';
    scen_cost.form='full';
    scen_cost.uels = {subs.uels};
    scen_output = rgdx('results',scen_cost);
    problem.optimal_scen_op_cost = scen_output.val;
    
    
    scen_run_idx = (1:scen_n)+(scen_n*(d_idx-1));
    problem.scen_op_cost(scen_run_idx) = scen_output.val;
    
end
end














