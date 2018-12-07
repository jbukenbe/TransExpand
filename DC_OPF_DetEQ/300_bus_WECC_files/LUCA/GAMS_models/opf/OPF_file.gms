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
x(gens, subs)                              Active power dispatched by generator g [MW]
pns(i, subs)                            Non-served power in node i [MW]
sg(gens, subs)                             Generation shed at generator g [MW]
*loss(lines,i,j, subs)                   Losses in the line i-j [MW]
;

Equations
eObFun                                  Objective Function
eSubFun(subs)                           Subset Objective Function
eDemand(i, subs)                        Demand met from generation and power inflow
ePowerFlow(lines,i,j, subs)             DC power flow in existing lines
*eLineLoss(lines,i,j,k,subs)             Line loss for existing lines
;

eObFun..                               z =e= sum(subs, zs(subs))/sum(subs,1);

eSubFun(subs)..                        zs(subs) =e= sum(gens, pVarCost(gens,subs) * x(gens, subs) + csg * sg(gens, subs))
                                           + sum(i, pns(i, subs) * cpns);


eDemand(i, subs)..                      sum(gens$ig(i,gens), x(gens, subs) - sg(gens, subs))
                                            - sum((lines,j) $ ii(lines,i,j), p(lines,i,j,subs))
                                            + sum((lines,j) $ ii(lines,j,i), p(lines,j,i, subs))
                                            =e= pLoad(i, subs) - pns(i, subs);

*eDemand(i, subs)..                      sum(g$ig(i,g), x(g, subs) - sg(g,subs))
*                                            - sum((lines,j) $ ii(lines,i,j), p(lines,i,j,subs)+loss(lines,i,j,subs)/2)
*                                            + sum((lines,j) $ ii(lines,j,i), p(lines,j,i, subs)-loss(lines,j,i,subs)/2)
*                                            =g= pLoad(i, subs) - pns(i, subs);


ePowerFlow(lines,i,j,subs)$ii(lines,i,j)..             p(lines,i,j, subs) =e= 1000*(theta(i, subs) - theta(j, subs)) * (dtl(lines,i,j,'s'));

*eLineLoss(lines,i,j,k,subs)$ii(lines,i,j)..            loss(lines,i,j,subs) =g= ((2*dtl(lines,i,j,'r'))/(dtl(lines,i,j,'r')**2 + (1/dtl(lines,i,j,'s'))**2))*((1000/300)*(theta(i,subs) - theta(j,subs)) * ds(k,'m') + ds(k,'n'));


ii(lines,i,j) $ dtl(lines,i,j,'pmax') = YES;
x.up(gens, subs) = pgen_max(gens);
x.lo(gens, subs) = pgen_min(gens);
p.up(lines,i,j, subs)$ii(lines,i,j) = dtl(lines,i,j,'pmax');
p.lo(lines,i,j, subs)$ii(lines,i,j) = -dtl(lines,i,j,'pmax');
theta.up(i, subs) = thetamax;
theta.lo(i, subs) = -thetamax;

x.fx(gens, subs)$(nd(gens)) = pRenew(gens,subs)*pgen_max(gens);
theta.fx(i, subs) $ [ORD(i) = 1] = 0;

Model opf /all/;

option optCR = .01;
option solprint = off;
option limrow = 0;
option limcol = 0;
option threads = 4;
option LP = Cplex;

opf.OptFile = 1

*option dmpsym;
Solve opf using LP minimizing z;

*$include "OpfOut.INC"


Parameters
scen_cost(subs)        cost of each hour
p_shed_gen(g,subs)  shed generation at each hour
p_pns(i,subs)          power not served at each hour
;
scen_cost(subs) = zs.l(subs);
p_shed_gen(gens,subs) = sg.l(gens, subs);
p_pns(i,subs) = pns.l(i,subs);


execute_unload 'results', scen_cost, p_shed_gen, p_pns;
