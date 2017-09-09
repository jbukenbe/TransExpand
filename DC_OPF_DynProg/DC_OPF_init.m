function params = DC_OPF_init(params)
% This file initializes the parameters for the DC optimal power flow operations
% model. It currently initializes the IEEE 30-Bus test case

%History            
%Version    Date        Who     Summary
%1          09/09/2017  JesseB  Initial Version

            
%% Scalor Data
params.theta_lim = .5;
params.cpns = 1000;
params.fix_line_cost = 4;
params.var_line_cost = 100;

%% Scenario Initialization
params.scen.n = 9;

%% Bus Initialization
params.bus.list = bus_list;
params.bus.load = load_list;
params.bus.n = size(params.bus.list,1);


%% Generator Initialization
params.gen.list
params.gen.loc
params.gen.p_max
params.gen.p_min
params.gen.op_cost
params.gen.n = size(params.gen.list,1);

%% Line Initialization
params.line.built
params.line.cand
params.line.curent
params.line.full
params.line.res
params.line.imp
params.line.max_flow
params.line.adj
params.line.from_to


params.line.n = size(params.line.full,1);

%% Line Loss Estimator Initialization
params.line.loss.m
params.line.loss.k
params.line.loss.size = size(params.line.loss.m,1);

end
