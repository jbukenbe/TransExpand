$ontext
This model solves the transmission expansion plan that minimizes the cost of the
DC OPF model for the WECC 300-Bus test case.

Options:
         Active Scenarios: Choose which scenarios to include in analysis e.g. /s2/ or /s1*s3/ etc...


History
Version__Date____________Who_____Summary
1        03/06/2018      JesseB  Initial Model adapted from IEEE 30 bus case


TODO:
         Add renewable limits

$offtext


$include        "WECC_Gams.txt"

SET subs(s)      Active Scenarios       ;

$set matout "'results.gdx',scen_cost";
$if not set gdxin $set gdxin MtoG
$GDXIN %gdxin%
$LOAD  subs
$GDXIN



SCALARS
ThetaMax         Maximum value of theta [radians] /.5/
cpns             Cost of non-served power         /1000/
csg              Cost of shed renewable gen       /300/
;


VARIABLES
z                Value of the objective function
p(l, i, j, subs)           Active power Flow in line i - j [MW]
theta(i, subs)         Angle of voltage in node i [rad]
;

POSITIVE VARIABLES
x(g, subs)             Active power dispatched by generator g [MW]
pns(i, subs)           Non-served power in node i [MW]
sg(g, subs)            Generation shed at generator g [MW]
loss(l,i,j, subs)      Losses in the line i-j [MW]
;

Binary VARIABLES
y(l,i,j)            1 if line is built from i to j 0 otherwise
;

Equations
eSubFun                          Objective Function
eDemand(i, subs)                    Demand met from generation and power inflow
eRenewGen(nd,subs)                  Renewable Generation limit
ePowerFlow(l,i,j, subs)             DC power flow in existing lines
eCandLineFlowHi(l,i,j,subs)         DC power flow upper bound in candidate lines
eCandLineFlowLow(l,i,j,subs)        DC power flow lower bound in candidate lines
eCandCapHi(l,i,j,subs)              Line capacity upper bound in candidate lines
eCandCapLow(l,i,j,subs)             Line capacity lower bound in candidate lines
eLineLoss(l,i,j,k,subs)             Line loss for existing lines
eCandLineLossHi(l,i,j,k,subs)       Line loss for candidate lines upper bound
eCandLineLossLow(l,i,j,k,subs)      Line loss for candidate lines lower bound
*etest                               line number test
;

eSubFun..                                           z =e= sum(subs, sum(g, pVarCost(g, subs) * x(g, subs) + csg * sg(g, subs))
                                                          + sum(i, pns(i, subs) * cpns))/sum(subs,1)
                                                          + sum((l,i,j)$cand(l,i,j),y(l,i,j)*dtl(l,i,j,'cost')*10);

eDemand(i, subs)..                                  sum(g$ig(i,g), x(g, subs) - sg(g,subs))
                                                          - sum((l,j) $ ii(l,i,j), p(l,i,j,subs)+loss(l,i,j,subs)/2)
                                                          + sum((l,j) $ ii(l,j,i), p(l,j,i, subs)-loss(l,j,i,subs)/2)
                                                          =g= pLoad(i, subs) - pns(i, subs);

eRenewGen(nd,subs)..                                x(nd, subs) =l= pRenew(nd,subs);

ePowerFlow(l,i,j,subs)$(built(l,i,j))..             p(l,i,j, subs) =e= 1000*(theta(i, subs) - theta(j, subs)) / (dtl(l,i,j,'x'));
eCandLineFlowHi(l,i,j,subs)$(cand(l,i,j))..         p(l,i,j, subs) =l= 10000*(1-y(l,i,j))+ 1000*(theta(i, subs) - theta(j, subs)) / (dtl(l,i,j,'x'));
eCandLineFlowLow(l,i,j,subs)$(cand(l,i,j))..        p(l,i,j, subs) =g= -10000*(1-y(l,i,j))+ 1000*(theta(i, subs) - theta(j, subs)) / (dtl(l,i,j,'x'));

eCandCapHi(l,i,j,subs)$(cand(l,i,j))..              p(l,i,j,subs) =l= y(l,i,j)*dtl(l,i,j,'pmax');
eCandCapLow(l,i,j,subs)$(cand(l,i,j))..             p(l,i,j,subs) =g= -y(l,i,j)*dtl(l,i,j,'pmax');

eLineLoss(l,i,j,k,subs)$(built(l,i,j))..            loss(l,i,j,subs) =g= ((2*dtl(l,i,j,'r'))/(dtl(l,i,j,'r')**2 + dtl(l,i,j,'x')**2))*((1000/300)*(theta(i,subs) - theta(j,subs)) * ds(k,'m') + ds(k,'n'));
eCandLineLossHi(l,i,j,k,subs)$(cand(l,i,j))..       loss(l,i,j,subs) =l= 1000*y(l,i,j);
eCandLineLossLow(l,i,j,k,subs)$(cand(l,i,j))..      loss(l,i,j,subs) =g= -1000*(1-y(l,i,j)) + ((2*dtl(l,i,j,'r'))/(dtl(l,i,j,'r')**2 + dtl(l,i,j,'x')**2))*((1000/300)*(theta(i,subs) - theta(j,subs)) * ds(k,'m') + ds(k,'n'));
*etest..                                             sum((l,i,j)$cand(l,i,j), y(l,i,j)) =l= 30;


ii(l,i,j) $ dtl(l,i,j,'pmax') = YES;
built(l,i,j)$dtl(l,i,j,'current')=YES;
cand(l,i,j) $dtl(l,i,j,'can')= YES;
x.up(g, subs) = pMax(g, subs);
x.lo(g, subs) = pMin(g, subs);
p.up(l,i,j, subs) = dtl(l,i,j,'pmax');
p.lo(l,i,j, subs) = -dtl(l,i,j,'pmax');
theta.up(i, subs) = thetamax;
theta.lo(i, subs) = -thetamax;


theta.fx(i, subs) $ [ORD(i) = 1] = 0;

*Model Master_DC_OPF /eMasterFun/
*Model Sub_DC_OPF /eSubFun, eDemand, ePowerFlow, ePCandFlowHi, ePCandFlowLo, ePCandLimitHi, ePCandLimitLo, ePGenCap, ePGenMin/;
Model opf /all/;
option optCR = .051;
option ResLim = 2000;
Solve opf using MIP minimizing z;

*$include "OpfOut.INC"

Parameter scen_cost(subs);
scen_cost(subs) = sum((g),x.l(g,subs));
execute_unload %matout%;