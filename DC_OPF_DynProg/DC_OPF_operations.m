function op_cost   = DC_OPF_operations(params, dec_var)
% This function optimizes the DC optimal power flow problem with line
% losses using the linprog function

%History            
%Version    Date        Who     Summary
%1          09/09/2017  JesseB  Initial version updated from old files


%% Data Initialization
theta_max = params.theta_lim;
line_n = sum(dec_var) + params.line.n;
bus_n = params.bus.n;
gen_n = params.gen.n;
m = 2*line_n;
n = gen_n + 2*bus_n + line_n;
A = zeros(m, n); 

%% Initialization of linear program structure
% Create operating cost segment
c_op = params.gen.op_cost;
% Create power not served cost segment
c_pns = params.cpns(ones(bus_n,1));
% Merge into full cost vector
c_free = zeros(1, (n-gen_n-bus_n));
c = [c_op, c_pns, c_free];


% Load constraints
% power not served matrix segment
A_pns = eye(bus_n);

% generation matrix segment
A_gen = zeros(bus_n, gen_n);
for g_idx = 1:gen_n
    A_gen(params.gen.loc(g_idx),g_idx) = 1;
end

% power flow matrix segment
A_flow = zeros(bus_n)
for l_idx = 1:line_n
    A_flow(params.line.from_to(l_idx,1),l_idx
end

% line loss matrix segment
A_loss = zeros(bus_n, line_n);
for l_idx = 1:line_n
    A_loss(params.line.from_to(l_idx,1),l_idx)= -.5;
    A_loss(params.line.from_to(l_idx,2),l_idx)= -.5;
end

% b vector
b_demand = params.bus.load;





% DC Power flow Constraints

% Line Loss Constraints


% Variable upper and lower bounds



% Transmission Capacity Constraints

PF = zeros(line_n,bus_n);
for i = 1:line_n
    PF(i,T(i,1)) = 1/params.line.imp(i);
    PF(i,T(i,2)) = -1/X(i);
end      
A(1:m,gen_n+1:gen_n+bus_n) = cat(1,PF,-PF);
b = repmat(pmax,2,1);


% The system Aeq*x = beq only covers demand constraints
DRec = zeros(bus);
for i = 1:line
    DRec(T(i,1),T(i,2)) =  -133/(X(i));
    DRec(T(i,2),T(i,1)) =  -133/(X(i));
end
DRec = DRec - diag(sum(DRec,2));
Aeq(1:bus,gen_n+1:gen_n+bus) = DRec;
Aeq(1:bus,gen_n+bus+1:gen_n+2*bus) = eye(bus);
beq = L;

% The rest is covered by variable bounds
ub = zeros(1,n); lb = zeros(1,n);
ub(1:gen_n) = GUp; ub(gen_n+1:gen_n+bus) = theta_max; ub(gen_n+bus+1:end) = inf;
lb(1:gen_n) = GLo; lb(gen_n+1:gen_n+bus) = -theta_max;



%% Solve the program
[~, op_cost] = linprog(c, A, b, Aeq, beq, lb, ub);

end

