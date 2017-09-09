function y = DC_OPF(problem)
% This function optimizes the DC optimal power flow problem


ThetaMax = problem.params.theta_lim;

line = size(T,1); bus = size(L,1); gen = size(GNode,1);
m = 2*line; n = gen + 2*bus;

c = zeros(1,n); c(1:gen) = OpCost; c(gen+bus+1:end) = 1000;

% The system A*x <= b only covers transmission capacity constraints
A = zeros(m, n); 
PF = zeros(line,bus);
for i = 1:line
    PF(i,T(i,1)) = 1/X(i); PF(i,T(i,2)) = -1/X(i);
end      
A(1:m,gen+1:gen+bus) = cat(1,PF,-PF); b = repmat(pmax,2,1);

% The system Aeq*x = beq only covers demand constraints
Aeq = zeros(bus,n);
for i = 1:gen
    Aeq(GNode(i),i) = 1;
end
DRec = zeros(bus);
for i = 1:line
    DRec(T(i,1),T(i,2)) =  -133/(X(i)); DRec(T(i,2),T(i,1)) =  -133/(X(i));
end
DRec = DRec - diag(sum(DRec,2));
Aeq(1:bus,gen+1:gen+bus) = DRec;
Aeq(1:bus,gen+bus+1:gen+2*bus) = eye(bus);
beq = L;

% The rest is covered by variable bounds
ub = zeros(1,n); lb = zeros(1,n);
ub(1:gen) = GUp; ub(gen+1:gen+bus) = ThetaMax; ub(gen+bus+1:end) = inf;
lb(1:gen) = GLo; lb(gen+1:gen+bus) = -ThetaMax;

[~, y] = linprog(c, A, b, Aeq, beq, lb, ub);

end

