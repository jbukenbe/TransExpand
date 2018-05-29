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
zh(t,subs,subh)               Value of the objective function in each possible hour
p(l, i, j, t, subs, subh)     Active power Flow in line i - j [MW]
theta(i, t, subs, subh)       Angle of voltage in node i [rad]
;

POSITIVE VARIABLES
x(g, t, subs, subh)    Active power dispatched by generator g [MW]
pns(i, t, subs, subh)  Non-served power in node i [MW]
sg(g, t, subs, subh)   Generation shed at generator g [MW]
*loss(l,i,j, subs)      Losses in the line i-j [MW]
;

Binary VARIABLES
y(cand,i,j,t,subs)            1 if line is built from i to j 0 otherwise
;

Equations
eObFun                                      Objective Function
eSubFun(t, subs, subh)                      Objective function at each scenario
eLineExist(cand,i,j,t,subs)                 Lines built must exist in future time stages CURRENTLY ONLY WORKS FOR 2 INVESMTENT STAGES
eQuickFix(cand,i,j,t,subs)                  First time stage has only one scenario
eDemand(i, t, subs, subh)                   Demand met from generation and power inflow
ePowerFlow(built,i,j,t, subs, subh)         DC power flow in existing lines
eCandLineFlowHi(cand,i,j,t,subs, subh)      DC power flow upper bound in candidate lines
eCandLineFlowLow(cand,i,j,t,subs, subh)     DC power flow lower bound in candidate lines
eCandCapHi(cand,i,j,t, subs, subh)          Line capacity upper bound in candidate lines
eCandCapLow(cand,i,j,t, subs, subh)         Line capacity lower bound in candidate lines
;

eObFun..                               z =e= sum(t,.5*sum(subs,psubs_w(subs)*(sum(subh, psubh_w(subh)*zh(t,subs,subh))
                                          + sum((cand,i,j)$ii(cand,i,j),y(cand,i,j,t,subs)*dtl(cand,i,j,'cost')*.1))));


eSubFun(t,subs,subh)..                 zh(t,subs,subh)
                                           =e= sum(g, pVarCost(g,t,subs,subh) * x(g, t,subs,subh) + csg * sg(g, t,subs,subh))
                                           + sum(i, pns(i, t,subs,subh) * cpns);


eLineExist(cand,i,j,t,subs)$((ORD(t)>1)
           and ii(cand,i,j))..          y(cand,i,j,t,subs) =g= y(cand,i,j,t-1,subs);


eQuickFix(cand,i,j,t,subs)$((ORD(t) = 1)
          and ii(cand,i,j))..           y(cand,i,j,t,subs) =e= y(cand,i,j,t,'s1');


eDemand(i, t,subs,subh)..              sum(g$ig(i,g), x(g,t,subs,subh) - sg(g,t,subs,subh))
                                           - sum((l,j) $ ii(l,i,j), p(l,i,j,t,subs,subh))
                                           + sum((l,j) $ ii(l,j,i), p(l,j,i, t,subs,subh))
                                           =e= pLoad(i, t,subs,subh) - pns(i, t,subs,subh);

*eDemand(i, subs)..                      sum(g$ig(i,g), x(g, subs) - sg(g,subs))
*                                            - sum((lines,j) $ ii(lines,i,j), p(lines,i,j,subs)+loss(lines,i,j,subs)/2)
*                                            + sum((lines,j) $ ii(lines,j,i), p(lines,j,i, subs)-loss(lines,j,i,subs)/2)
*                                            =g= pLoad(i, subs) - pns(i, subs);


ePowerFlow(built,i,j,t,subs,subh)$ii(built,i,j)..            p(built,i,j, t,subs,subh) =e= 1000*(theta(i, t,subs,subh) - theta(j, t,subs,subh)) * (dtl(built,i,j,'s'));
eCandLineFlowHi(cand,i,j,t,subs,subh)$ii(cand,i,j)..         p(cand,i,j, t,subs,subh) =l= 10000*(1-y(cand,i,j,t,subs))+ 1000*(theta(i, t,subs,subh) - theta(j, t,subs,subh)) * (dtl(cand,i,j,'s'));
eCandLineFlowLow(cand,i,j,t,subs,subh)$ii(cand,i,j)..        p(cand,i,j, t,subs,subh) =g= -10000*(1-y(cand,i,j,t,subs))+ 1000*(theta(i, t,subs,subh) - theta(j, t,subs,subh)) * (dtl(cand,i,j,'s'));

eCandCapHi(cand,i,j,t,subs,subh)$ii(cand,i,j)..              p(cand,i,j,t,subs,subh) =l= y(cand,i,j,t,subs)*dtl(cand,i,j,'pmax');
eCandCapLow(cand,i,j,t,subs,subh)$ii(cand,i,j)..             p(cand,i,j,t,subs,subh) =g= -y(cand,i,j,t,subs)*dtl(cand,i,j,'pmax');




ii(l,i,j) $ dtl(l,i,j,'pmax') = YES;
x.up(g, t,subs,subh) = gdata(g,'pmax');
x.lo(g, t,subs,subh) = gdata(g,'pmin');
p.up(l,i,j, t,subs,subh)$ii(l,i,j) = dtl(l,i,j,'pmax');
p.lo(l,i,j, t,subs,subh)$ii(l,i,j) = -dtl(l,i,j,'pmax');
theta.up(i, t,subs,subh) = thetamax;
theta.lo(i, t,subs,subh) = -thetamax;

x.fx(nd, t,subs,subh) = pRenew(nd,t,subs,subh)*gdata(nd,'pmax');
theta.fx(i, t,subs,subh)$([ORD(i) = 1]) = 0;

*Model Master_DC_OPF /eMasterFun/
*Model Sub_DC_OPF /eSubFun, eDemand, ePowerFlow, ePCandFlowHi, ePCandFlowLo, ePCandLimitHi, ePCandLimitLo, ePGenCap, ePGenMin/;
Model opf /all/;


* limit run time
option ResLim = 3000;

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
scen_cost(t,subs,subh)
tot_cost
lines_built(cand,i,j,t,subs)
;

scen_cost(t,subs,subh) = zh.l(t,subs,subh);
tot_cost = z.l;
lines_built(cand,i,j,t,subs)$ii(cand,i,j) = y.l(cand,i,j,t,subs);
execute_unload 'results', scen_cost, tot_cost, lines_built;
