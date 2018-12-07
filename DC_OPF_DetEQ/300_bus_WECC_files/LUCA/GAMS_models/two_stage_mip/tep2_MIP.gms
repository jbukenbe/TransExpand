$ontext
This model solves the transmission expansion plan that minimizes the cost of the
DC OPF model for the WECC 300-Bus test case.



History
Version__Date____________Who_____Summary
1        03/06/2018      JesseB  Initial Model adapted from IEEE 30 bus case
2        05/18/2018      JesseB  Adapted to take Matlab input
3        05/23/2018      JesseB  Expanded to 2 investment stages

TODO:
         Add renewable limits

$offtext
$include        "WECC_Gams_tep2.txt"

SCALARS
ThetaMax         Maximum value of theta [radians] /.5/
cpns             Cost of non-served power         /1000/
csg              Cost of shed renewable gen       /300/
r                Discount rate per time step       /.95/
;

VARIABLES
z                             Value of the objective function
zh(subs,subh)               Value of the objective function in each possible hour
p(l, i, j, subs, subh)     Active power Flow in line i - j [MW]
theta(i, subs, subh)       Angle of voltage in node i [rad]
;

POSITIVE VARIABLES
x(g, subs, subh)    Active power dispatched by generator g [MW]
pns(i, subs, subh)  Non-served power in node i [MW]
sg(g, subs, subh)   Generation shed at generator g [MW]
*loss(l,i,j, subs)      Losses in the line i-j [MW]
;

Binary VARIABLES
y(cand,i,j,subs)            1 if line is built from i to j 0 otherwise
;

Equations
eObFun                                      Objective Function
eSubFun(subs, subh)                      Objective function at each scenario
eLineExist(cand,i,j,subs)                 Lines built must exist in future time stages CURRENTLY ONLY WORKS FOR 2 INVESMTENT STAGES
eDemand(i, subs, subh)                   Demand met from generation and power inflow
ePowerFlow(built,i,j, subs, subh)         DC power flow in existing lines
eCandLineFlowHi(cand,i,j,subs, subh)      DC power flow upper bound in candidate lines
eCandLineFlowLow(cand,i,j,subs, subh)     DC power flow lower bound in candidate lines
eCandCapHi(cand,i,j, subs, subh)          Line capacity upper bound in candidate lines
eCandCapLow(cand,i,j, subs, subh)         Line capacity lower bound in candidate lines
;

eObFun..                               z =e= sum(subs,psubs_w(subs)*(sum(subh, psubh_w(subh)*zh(subs,subh))
                                          + sum((cand,i,j)$ii(cand,i,j),y(cand,i,j,subs)*dtl(cand,i,j,'cost'))));


eSubFun(subs,subh)..                 zh(subs,subh)
                                           =e= sum(g, pVarCost(g,subs,subh) * x(g,subs,subh) + csg * sg(g, subs,subh))
                                           + sum(i, pns(i, subs,subh) * cpns);


eLineExist(cand,i,j,subs)$(ii(cand,i,j)
            and ss(subs,subq))..            y(cand,i,j,subs) - y(cand,i,j,subq) =l= 0;


eDemand(i, subs,subh)..              sum(g$ig(i,g), x(g,subs,subh) - sg(g,subs,subh))
                                           - sum((l,j) $ ii(l,i,j), p(l,i,j,subs,subh))
                                           + sum((l,j) $ ii(l,j,i), p(l,j,i, subs,subh))
                                           =e= pLoad(i, subs,subh) - pns(i, subs,subh);


ePowerFlow(built,i,j,subs,subh)$ii(built,i,j)..            p(built,i,j, subs,subh) =e= 1000*(theta(i, subs,subh) - theta(j, subs,subh)) * (dtl(built,i,j,'s'));


eCandLineFlowHi(cand,i,j,subs,subh)$ii(cand,i,j)..         p(cand,i,j, subs,subh) =l= 10000*(1-y(cand,i,j,subs))+ 1000*(theta(i,subs,subh) - theta(j,subs,subh)) * (dtl(cand,i,j,'s'));
eCandLineFlowLow(cand,i,j,subs,subh)$ii(cand,i,j)..        p(cand,i,j, subs,subh) =g= -10000*(1-y(cand,i,j,subs))+ 1000*(theta(i,subs,subh) - theta(j,subs,subh)) * (dtl(cand,i,j,'s'));


eCandCapHi(cand,i,j,subs,subh)$ii(cand,i,j)..              p(cand,i,j,subs,subh) =l= y(cand,i,j,subs)*dtl(cand,i,j,'pmax');
eCandCapLow(cand,i,j,subs,subh)$ii(cand,i,j)..             p(cand,i,j,subs,subh) =g= -y(cand,i,j,subs)*dtl(cand,i,j,'pmax');




ii(l,i,j) $ dtl(l,i,j,'pmax') = YES;
x.up(g, subs,subh) = gdata(g,'pmax');
x.lo(g, subs,subh) = gdata(g,'pmin');
p.up(l,i,j, subs,subh)$ii(l,i,j) = dtl(l,i,j,'pmax');
p.lo(l,i,j, subs,subh)$ii(l,i,j) = -dtl(l,i,j,'pmax');
theta.up(i, subs,subh) = thetamax;
theta.lo(i, subs,subh) = -thetamax;

x.fx(nd, subs,subh) = pRenew(nd,subs,subh)*gdata(nd,'pmax');
theta.fx(i, subs,subh)$([ORD(i) = 1]) = 0;

*Model Master_DC_OPF /eMasterFun/
*Model Sub_DC_OPF /eSubFun, eDemand, ePowerFlow, ePCandFlowHi, ePCandFlowLo, ePCandLimitHi, ePCandLimitLo, ePGenCap, ePGenMin/;
Model opf /all/;


* limit run time
option ResLim = 1500;

* optimality gap
option optCR = .12;

* print solution file
option solprint = on;

* shorten output colums and rows
option limrow = 0;
option limcol = 0;

* run in parallel with  threads
option threads = 4;

* use Cplex solver for linear programs
option LP = Cplex;

* use OPT file cplex for minimizing lp problems concurrently
opf.OptFile = 1

option dmpsym;
Solve opf using MIP minimizing z;

*$include "OpfOut.INC"



Parameters
scen_cost(subs,subh)
tot_cost
lines_built(cand,i,j,subs)
;

scen_cost(subs,subh) = zh.l(subs,subh);
tot_cost = z.l;
lines_built(cand,i,j,subs)$ii(cand,i,j) = y.l(cand,i,j,subs);
execute_unload 'results', scen_cost, tot_cost, lines_built;
