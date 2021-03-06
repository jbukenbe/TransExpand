$ontext
This model optimizes the DC OPF problem for use in transmission expansion planning
The lines used in this model must be specified by a matlab interface


History
Version__Date____________Who_____Summary
1        03/10/2018      JesseB  Initial Model adapted Transmission Expansion model


$offtext

$offlisting
$include        "WECC_Gams.txt"
*$onlisting


SCALARS
ThetaMax         Maximum value of theta [radians] /.5/
cpns             Cost of non-served power         /1000/
csg              Cost of shed renewable gen       /300/
;


VARIABLES
z                          Value of the objective function
zs(subs)                   value of the objective function in each scenario
p(lines, i, j, subs)       Active power Flow in line i - j [MW]
theta(i, subs)             Angle of voltage in node i [rad]
;

POSITIVE VARIABLES
x(g, subs)                              Active power dispatched by generator g [MW]
pns(i, subs)                            Non-served power in node i [MW]
sg(g, subs)                             Generation shed at generator g [MW]
loss(lines,i,j, subs)                   Losses in the line i-j [MW]
;

Equations
eObFun                                  Objective Function
eSubFun(subs)                           Subset Objective Function
eDemand(i, subs)                        Demand met from generation and power inflow
eRenewGen(nd,subs)                      Renewable Generation limit
ePowerFlow(lines,i,j, subs)             DC power flow in existing lines
eLineLoss(lines,i,j,k,subs)             Line loss for existing lines
;

eObFun..                               z =e= sum(subs, zs(subs))/sum(subs,1);

eSubFun(subs)..                        zs(subs) =e= sum(g, pVarCost(g,subs) * x(g, subs) + csg * sg(g, subs))
                                           + sum(i, pns(i, subs) * cpns);

eDemand(i, subs)..                      sum(g$ig(i,g), x(g, subs) - sg(g,subs))
                                            - sum((lines,j) $ ii(lines,i,j), p(lines,i,j,subs)+loss(lines,i,j,subs)/2)
                                            + sum((lines,j) $ ii(lines,j,i), p(lines,j,i, subs)-loss(lines,j,i,subs)/2)
                                            =g= pLoad(i, subs) - pns(i, subs);

eRenewGen(nd,subs)..                     x(nd, subs) =e= pRenew(nd,subs)*gdata(nd,'pmax');

ePowerFlow(lines,i,j,subs)$ii(lines,i,j)..             p(lines,i,j, subs) =e= 1000*(theta(i, subs) - theta(j, subs)) * (dtl(lines,i,j,'s'));

eLineLoss(lines,i,j,k,subs)$ii(lines,i,j)..            loss(lines,i,j,subs) =g= ((2*dtl(lines,i,j,'r'))/(dtl(lines,i,j,'r')**2 + (1/dtl(lines,i,j,'s'))**2))*((1000/300)*(theta(i,subs) - theta(j,subs)) * ds(k,'m') + ds(k,'n'));


ii(lines,i,j) $ dtl(lines,i,j,'pmax') = YES;
x.up(g, subs) = gdata(g,'pmax');
x.lo(g, subs) = gdata(g,'pmin');
p.up(lines,i,j, subs)$ii(lines,i,j) = dtl(lines,i,j,'pmax');
p.lo(lines,i,j, subs)$ii(lines,i,j) = -dtl(lines,i,j,'pmax');
theta.up(i, subs) = thetamax;
theta.lo(i, subs) = -thetamax;


theta.fx(i, subs) $ [ORD(i) = 1] = 0;

Model opf /all/;
option optCR = .051;
option solprint = off;
option limrow = 0;
option limcol = 0;
*option dmpsym;
Solve opf using LP minimizing z;

*$include "OpfOut.INC"

Parameter scen_cost(subs);
scen_cost(subs) = zs.l(subs);
Display scen_cost;
execute_unload 'results', scen_cost;
