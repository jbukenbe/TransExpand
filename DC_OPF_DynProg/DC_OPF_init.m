function params = DC_OPF_init()
% This file initializes the parameters for the DC optimal power flow operations
% model. It currently initializes the IEEE 30-Bus test case

%History            
%Version    Date        Who     Summary
%1          09/09/2017  JesseB  Initial Version
%2          10/07/2017  JesseB  Added PLS regression parameters
%3          10/08/2017  JesseB  Added line cost data

%% Read GAMS Data
bus_data = importdata('bus_data.txt');
gen_loc_data = importdata('gen_loc_data.txt');
gen_price_data = importdata('gen_price_data.txt');
gen_pmin_data = importdata('gen_pmin_data.txt');
gen_pmax_data = importdata('gen_pmax_data.txt');
line_data = importdata('line_data.txt');
cos_apx_data = importdata('cos_apx_data.txt');


%% PLS Data
params.pls.interaction = 1;
params.pls.refine_samp_n = 400;
params.pls.line_samp_n = 5000;
params.pls.fit_samp_n = 1000000;
params.pls.n_comp = 10;

%% Scalor Optimization Data
params.theta_lim = .5;
params.cpns = 1000;
params.fix_line_cost = 4;
params.var_line_cost = 100;
params.initial_samp_n = 500;

%% Scenario Initialization
params.scen.n = size(bus_data.data,2);
params.scen.p = ones(params.scen.n,1)*1/params.scen.n;

%% Bus Initialization
params.bus.load = bus_data.data;
params.bus.n = size(params.bus.load,1);


%% Generator Initialization
params.gen.loc = gen_loc_data.data;
params.gen.p_max = gen_pmax_data.data;
params.gen.p_min = gen_pmin_data.data;
params.gen.op_cost = gen_price_data.data;
params.gen.n = size(params.gen.loc,1);

%% Line Initialization
params.line.cost = line_data.data(:,8);
params.line.built = line_data.data(:,7);
params.line.cand = line_data.data(:,6);
params.line.full = line_data.data(:,6) + line_data.data(:,7);
params.line.res = line_data.data(:,3);
params.line.imp = line_data.data(:,4);
params.line.max_flow = line_data.data(:,5);
params.line.from_to = [line_data.data(:,1),line_data.data(:,2)];
params.line.n = size(params.line.full,1);

%% Line Loss Estimator Initialization
params.line.loss.slope = cos_apx_data.data(:,1);
params.line.loss.const = cos_apx_data.data(:,2);
params.line.loss.n = size(params.line.loss.slope,1);

end
