$ontext
This model solves the transmission expansion plan that minimizes the cost of the
DC OPF model for the WECC 300-Bus test case.

Options:
         Active Scenarios: Choose which scenarios to include in analysis e.g. /s2/ or /s1*s3/ etc...


History
Version__Date____________Who_____Summary
1        03/06/2018      JesseB  Initial Model adapted from IEEE 30 bus case
2        05/18/2018      JesseB  Adapted to take Matlab input

TODO:
         Add renewable limits

$offtext
$include        "WECC_Gams.txt"

SCALARS
ThetaMax         Maximum value of theta [radians] /.5/
cpns             Cost of non-served power         /1000/
csg              Cost of shed renewable gen       /300/
;

VARIABLES
z                     Value of the objective function
zs(subs)              value of the objective function in each scenario
p(l, i, j, subs)      Active power Flow in line i - j [MW]
theta(i, subs)        Angle of voltage in node i [rad]
;

POSITIVE VARIABLES
x(g, subs)             Active power dispatched by generator g [MW]
pns(i, subs)           Non-served power in node i [MW]
sg(g, subs)            Generation shed at generator g [MW]
*loss(l,i,j, subs)      Losses in the line i-j [MW]
;

Binary VARIABLES
y(cand,i,j)            1 if line is built from i to j 0 otherwise
;

Equations
eObFun                          Objective Function
eSubFun(subs)                   Objective function at each scenario
eDemand(i, subs)                    Demand met from generation and power inflow
ePowerFlow(built,i,j, subs)             DC power flow in existing lines
eCandLineFlowHi(cand,i,j,subs)         DC power flow upper bound in candidate lines
eCandLineFlowLow(cand,i,j,subs)        DC power flow lower bound in candidate lines
eCandCapHi(cand,i,j,subs)              Line capacity upper bound in candidate lines
eCandCapLow(cand,i,j,subs)             Line capacity lower bound in candidate lines
;

eObFun..                               z =e= sum(subs, psub_w(subs)*zs(subs))+ sum((cand,i,j)$ii(cand,i,j),y(cand,i,j)*dtl(cand,i,j,'cost'));

eSubFun(subs)..                        zs(subs) =e= sum(g, pVarCost(g,subs) * x(g, subs) + csg * sg(g, subs))
                                           + sum(i, pns(i, subs) * cpns);

eDemand(i, subs)..                     sum(g$ig(i,g), x(g, subs) - sg(g,subs))
                                           - sum((l,j) $ ii(l,i,j), p(l,i,j,subs))
                                           + sum((l,j) $ ii(l,j,i), p(l,j,i, subs))
                                           =e= pLoad(i, subs) - pns(i, subs);

*eDemand(i, subs)..                      sum(g$ig(i,g), x(g, subs) - sg(g,subs))
*                                            - sum((lines,j) $ ii(lines,i,j), p(lines,i,j,subs)+loss(lines,i,j,subs)/2)
*                                            + sum((lines,j) $ ii(lines,j,i), p(lines,j,i, subs)-loss(lines,j,i,subs)/2)
*                                            =g= pLoad(i, subs) - pns(i, subs);


ePowerFlow(built,i,j,subs)$ii(built,i,j)..            p(built,i,j, subs) =e= 1000*(theta(i, subs) - theta(j, subs)) * (dtl(built,i,j,'s'));
eCandLineFlowHi(cand,i,j,subs)$ii(cand,i,j)..         p(cand,i,j, subs) =l= 10000*(1-y(cand,i,j))+ 1000*(theta(i, subs) - theta(j, subs)) * (dtl(cand,i,j,'s'));
eCandLineFlowLow(cand,i,j,subs)$ii(cand,i,j)..        p(cand,i,j, subs) =g= -10000*(1-y(cand,i,j))+ 1000*(theta(i, subs) - theta(j, subs)) * (dtl(cand,i,j,'s'));

eCandCapHi(cand,i,j,subs)$ii(cand,i,j)..              p(cand,i,j,subs) =l= y(cand,i,j)*dtl(cand,i,j,'pmax');
eCandCapLow(cand,i,j,subs)$ii(cand,i,j)..             p(cand,i,j,subs) =g= -y(cand,i,j)*dtl(cand,i,j,'pmax');




ii(l,i,j) $ dtl(l,i,j,'pmax') = YES;
*built(l,i,j)$dtl(l,i,j,'current')=YES;
*cand(l,i,j) $dtl(l,i,j,'can')= YES;
x.up(g, subs) = gdata(g,'pmax');
x.lo(g, subs) = gdata(g,'pmin');
p.up(l,i,j, subs)$ii(l,i,j) = dtl(l,i,j,'pmax');
p.lo(l,i,j, subs)$ii(l,i,j) = -dtl(l,i,j,'pmax');
theta.up(i, subs) = thetamax;
theta.lo(i, subs) = -thetamax;

x.fx(nd, subs) = pRenew(nd,subs)*gdata(nd,'pmax');
theta.fx(i, subs) $ [ORD(i) = 1] = 0;

*Model Master_DC_OPF /eMasterFun/
*Model Sub_DC_OPF /eSubFun, eDemand, ePowerFlow, ePCandFlowHi, ePCandFlowLo, ePCandLimitHi, ePCandLimitLo, ePGenCap, ePGenMin/;
Model opf /all/;


* limit run time
option ResLim = 500;

* optimality gap
option optCR = .05;

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
scen_cost(subs)
tot_cost
lines_built(cand,i,j)
;

scen_cost(subs) = zs.l(subs);
tot_cost = z.l;
lines_built(cand,i,j)$ii(cand,i,j) = y.l(cand,i,j);
execute_unload 'results', scen_cost, tot_cost, lines_built;
