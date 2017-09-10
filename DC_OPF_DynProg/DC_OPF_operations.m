function op_cost   = DC_OPF_operations(params, dec_var)
% This function optimizes the DC optimal power flow problem with line
% losses using the linprog function

%History            
%Version    Date        Who     Summary
%1          09/09/2017  JesseB  Initial version updated from old files


%% Data Initialization
theta_lim = params.theta_lim;
line_n = sum(dec_var) + params.line.n;
bus_n = params.bus.n;
gen_n = params.gen.n;
cos_n = params.line.loss.n;


%% Cost Vector
% Create operating cost segment
c_gen = params.gen.op_cost';
% Create power not served cost segment
c_pns = params.cpns(ones(bus_n,1))';
% Merge into full cost vector
c_theta = zeros(1,bus_n);
c_loss = zeros(1,line_n);
c = [c_gen, c_pns, c_theta, c_loss];


%% Load constraints
% generation matrix segment
A_load_gen = zeros(bus_n, gen_n);
for g_idx = 1:gen_n
    A_load_gen(params.gen.loc(g_idx),g_idx) = -1;
end

% power not served matrix segment
A_load_pns = -eye(bus_n);

% power flow (theta) matrix segment
A_load_theta = zeros(bus_n);
for l_idx = 1:line_n
    i = params.line.from_to(l_idx,1);
    j = params.line.from_to(l_idx,2);
    A_load_theta(i,i) = A_load_theta(i,i) - 1/params.line.imp(l_idx);
    A_load_theta(i,j) = A_load_theta(i,j) + 1/params.line.imp(l_idx);
    A_load_theta(j,j) = A_load_theta(j,j) - 1/params.line.imp(l_idx);
    A_load_theta(j,i) = A_load_theta(j,i) + 1/params.line.imp(l_idx);
end

% line loss matrix segment
A_load_loss = zeros(bus_n, line_n);
for l_idx = 1:line_n
    A_load_loss(params.line.from_to(l_idx,1),l_idx)= .5;
    A_load_loss(params.line.from_to(l_idx,2),l_idx)= .5;
end
% b vector
b_load= -params.bus.load;


%% Line Loss Constraints
%generation and pns segments are all zeros
A_loss_gen = zeros(line_n*cos_n, gen_n);
A_loss_pns = zeros(line_n*cos_n, bus_n);

%cosine appriximation coefficients
A_loss_theta = zeros(line_n*cos_n, bus_n);
b_loss = zeros(line_n*cos_n,1); 
for l_idx =1:l_idx
    res = params.line.res;
    imp = params.line.imp;
    i = params.line.from_to(l_idx,1);
    j = params.line.from_to(l_idx,2);
    for k_idx = 1:cos_n
        const = params.line.loss.const(k_idx);
        slope = params.line.loss.slope(k_idx);
        loss_coef = slope(k_idx)*(imp^2)*(1000/133);
        a_idx = (l_idx-1)*cos_n+(k_idx-1)+1;
        A_loss_theta(a_idx, i)=  loss_coef;
        A_loss_theta(a_idx, j)=  -loss_coef;
        b_loss(a_idx) = -const+(2*res/(res^2 +imp^2));
    end
end

%loss variable coefficients
A_loss_loss = -eye(line_n);
A_loss_loss = repelem(A_loss_loss,1,cos_n);


%% Variable upper and lower bounds
ub_gen = params.gen.pmax;
ub_pns = ones(1,bus_n)*inf;
ub_theta = ones(1,bus_n)*theta_lim;
ub_loss = ones(1,line_n)*inf;
lb_gen = params.gen.pmin;
lb_pns = zeros(1,bus_n);
lb_theta = -ones(1,bus_n)*theta_lim;
lb_loss = zeros(1,line_n);

ub = [ub_gen, ub_pns, ub_theta, ub_loss];
lb = [lb_gen, lb_pns, lb_theta, lb_loss];


%% Problem Solution Phase
A = [A_load_gen, A_load_pns, A_load_theta, A_load_loss;...
    A_loss_gen, A_loss_pns, A_loss_theta, A_loss_loss];
b = [b_load; b_loss];

[~, op_cost] = linprog(c, A, b, lb, ub);

end

