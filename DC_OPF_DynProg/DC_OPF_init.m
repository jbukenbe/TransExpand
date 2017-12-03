function params = DC_OPF_init(problem)
% This file initializes the parameters for the DC optimal power flow operations
% model. It currently initializes the IEEE 30-Bus test case

%History            
%Version    Date        Who     Summary
%1          09/09/2017  JesseB  Initial Version
%2          10/07/2017  JesseB  Added PLS regression parameters
%3          10/08/2017  JesseB  Added line cost data
%4          11/18/2017  JesseB  Added Parameters for online search
%5          11/30/2017  JesseB  Added maximim number of new lines to install
%6          12/02/2017  JesseB  Added candidate line info

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
params.pls.line_samp_n = 5000000;
params.pls.fit_samp_n = 5000000;
params.pls.n_comp = 10;
params.pls.dist_mat_size = 100000;

%% Scalor Optimization Data
params.theta_lim = .5;
params.cpns = 1000;
params.fix_line_cost = 4;
params.var_line_cost = 100;
params.max_new_lines = 23;
params.line_budget = 150;
params.initial_samp_n = 200;
params.refine_samp_n = 130;

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

%% Load Candidate Line Info

% pick candidate line ids
switch problem.problem_size
    case 196
        line_id = (42:237)';
    case 75
        line_id = [42 43 50 51 52 53 55 57 60 62 63 64 65 68 69 70 71 72 73 75 76 78 80 82 84 85 86 92 94 95 96 98 100 101 103 105 108 110 112 114 116 118 120 122 125 126 128 130 131 140:164 237]';
    case 50
        line_id = [42 43 50 51 52 53 55 57 60 62 63 64 65 68 69 70 71 72 73 75 76 78 80 82 84 85 86 92 94 95 96 98 100 101 103 105 108 110 112 114 116 118 120 122 125 126 128 130 131 237]';
    case 45
        line_id = [42 50 51 52 53 55 57 60 62 63 64 65 68 69 71 72 73 75 76 78 80 82 84 85 86 92 94 96 98 100 101 103 105 108 112 114 116 118 120 122 125 126 128 130 131]';
    case 40
        line_id = [42 50 51 52 53 55 57 60 62 63 64 65 68 69 72 73 75 78 80 82 84 85 86 92 94 96 98 100 101 103 105 108 112 114 116 118 120 122 125 131]';
    case 35
        line_id = [50 51 52 53 55 57 60 62 63 64 65 68 69 72 73 75 78 80 82 84 85 86 92 94 96 98 100 101 103 105 108 112 118 120 131]';
    case 30
        line_id = [50 51 52 53 55 57 60 62 63 64 65 68 69 72 73 75 78 80 82 84 85 86 92 94 96 98 100 101 118 131]';
    case 25
        line_id = [51 52 53 55 57 60 62 63 64 65 68 69 72 75 82 84 85 86 92 94 96 98 100 118 131]';
    case 20
        line_id = [51 52 53 55 57 60 62 64 68 69 72 75 82 92 94 96 98 100 118 131]';
    case 15
        line_id = [51 52 53 55 57 60 62 64 69 72 75 82 94 118 131]';
    case 10
        line_id = [51 52 53 55 57 72 75 82 118 131]';
    otherwise
        line_id = (42:237)';
end

% write problem specific info based on included lines
params.cand.line_id{1} = line_id;
params.new_line_cost{1} = params.line.cost(line_id);
params.cand.n = length(line_id);
params.plan.n = 2^params.cand.n;
params.line.dec_built{1} = logical(params.line.built);
end
