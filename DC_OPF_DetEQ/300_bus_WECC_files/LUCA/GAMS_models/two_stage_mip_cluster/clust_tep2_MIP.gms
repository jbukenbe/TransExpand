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

*$offlisting
$include        "WECC_clust_tep2.txt"

SCALARS
ThetaMax         Maximum value of theta [radians] /.5/
cpns             Cost of non-served power         /1000/
csg              Cost of shed renewable gen       /300/
r                Discount rate per time step       /1/
;

VARIABLES
z                             Value of the objective function
zh(c,subs,subh)               Value of the objective function in each possible hour
p(l, i,j,c,subs,subh)     Active power Flow in line i - j [MW]
theta(i,c,subs,subh)       Angle of voltage in node i [rad]
;

POSITIVE VARIABLES
x(g,c,subs,subh)    Active power dispatched by generator g [MW]
pns(i,c,subs,subh)  Non-served power in node i [MW]
shg(g,c,subs,subh)   Generation shed at generator g [MW]
*loss(l,i,j, subs)      Losses in the line i-j [MW]
;

Binary VARIABLES
y(cand,i,j,c)            1 if line is built from i to j 0 otherwise
;

Equations
eObFun                                      Objective Function
eSubFun(c,subs,subh)                      Objective function at each scenario
eLineExist(cand,i,j,c,k)                 Lines built must exist in future time stages CURRENTLY ONLY WORKS FOR 2 INVESMTENT STAGES
eDemand(i,c,subs,subh)                   Demand met from generation and power inflow
ePowerFlow(built,i,j,c,subs,subh)         DC power flow in existing lines
eCandLineFlowHi(cand,i,j,c,subs,subh)      DC power flow upper bound in candidate lines
eCandLineFlowLow(cand,i,j,c,subs,subh)     DC power flow lower bound in candidate lines
eCandCapHi(cand,i,j,c,subs,subh)          Line capacity upper bound in candidate lines
eCandCapLow(cand,i,j,c,subs,subh)         Line capacity lower bound in candidate lines
;

*eObFun..                               z =e= sum(c,sum(s,psubs_w(s)*(sum(subh$csh(c,s,subh), psubh_w(subh)*zh(c,s,subh))
*                                          + sum((cand,i,j)$ii(cand,i,j),y(cand,i,j,c)*dtl(cand,i,j,'cost'))));


eObFun..                               z =e= sum(c, pc_w(c)*(pcluster_means(c) + sum((cand,i,j)$ii(cand,i,j),y(cand,i,j,c)*dtl(cand,i,j,'cost'))
                                          + sum((subs,subh)$csh(c,subs,subh), psub_w(c,subs,subh)*(zh(c,subs,subh)-ph_means(c,subs,subh)))));


eSubFun(c,subs,subh)$csh(c,subs,subh)..      zh(c,subs,subh)
                                           =e= sum(g$sg(subs,g), pVarCost(g,subs,subh) * x(g,c,subs,subh) + csg * shg(g,c,subs,subh))
                                           + sum(i, pns(i,c,subs,subh) * cpns);


eLineExist(cand,i,j,c,k)$(ii(cand,i,j)
            and cc(c,k))..             y(cand,i,j,c) - y(cand,i,j,k) =l= 0;


eDemand(i,c,subs,subh)$csh(c,subs,subh)..    sum(g$(ig(i,g) and sg(subs,g)), x(g,c,subs,subh) - shg(g,c,subs,subh))
                                           - sum((l,j) $ ii(l,i,j), p(l,i,j,c,subs,subh))
                                           + sum((l,j) $ ii(l,j,i), p(l,j,i,c,subs,subh))
                                           =e= pLoad(i,subs,subh) - pns(i,c,subs,subh);


ePowerFlow(built,i,j,c,subs,subh)$(ii(built,i,j)
           and csh(c,subs,subh))..        p(built,i,j,c,subs,subh) =e= 1000*(theta(i,c,subs,subh) - theta(j,c,subs,subh)) * (dtl(built,i,j,'s'));


eCandLineFlowHi(cand,i,j,c,subs,subh)$(ii(cand,i,j)
                and csh(c,subs,subh))..                          p(cand,i,j,c,subs,subh) =l= 10000*(1-y(cand,i,j,c))+ 1000*(theta(i,c,subs,subh) - theta(j,c,subs,subh)) * (dtl(cand,i,j,'s'));
eCandLineFlowLow(cand,i,j,c,subs,subh)$(ii(cand,i,j)
                and csh(c,subs,subh))..                          p(cand,i,j,c,subs,subh) =g= -10000*(1-y(cand,i,j,c))+ 1000*(theta(i,c,subs,subh) - theta(j,c,subs,subh)) * (dtl(cand,i,j,'s'));


eCandCapHi(cand,i,j,c,subs,subh)$(ii(cand,i,j)
                and csh(c,subs,subh))..                          p(cand,i,j,c,subs,subh) =l= y(cand,i,j,c)*dtl(cand,i,j,'pmax');
eCandCapLow(cand,i,j,c,subs,subh)$(ii(cand,i,j)
                and csh(c,subs,subh))..                          p(cand,i,j,c,subs,subh) =g= -y(cand,i,j,c)*dtl(cand,i,j,'pmax');




ii(l,i,j) $ dtl(l,i,j,'pmax') = YES;
x.up(g,c,subs,subh)$(sg(subs,g) and csh(c,subs,subh)) = pgen_max(g);
x.lo(g,c,subs,subh)$(sg(subs,g) and csh(c,subs,subh)) = pgen_min(g);
p.up(l,i,j,c,subs,subh)$(ii(l,i,j) and csh(c,subs,subh)) = dtl(l,i,j,'pmax');
p.lo(l,i,j,c,subs,subh)$(ii(l,i,j) and csh(c,subs,subh)) = -dtl(l,i,j,'pmax');
theta.up(i,c,subs,subh) = thetamax;
theta.lo(i,c,subs,subh) = -thetamax;

x.fx(nd,c,subs,subh)$(sg(subs,nd) and csh(c,subs,subh)) = pRenew(nd,subs,subh)*pgen_max(nd);
theta.fx(i,c,subs,subh)$(csh(c,subs,subh) and ([ORD(i) = 1])) = 0;

*Model Master_DC_OPF /eMasterFun/
*Model Sub_DC_OPF /eSubFun, eDemand, ePowerFlow, ePCandFlowHi, ePCandFlowLo, ePCandLimitHi, ePCandLimitLo, ePGenCap, ePGenMin/;



Model opf /all/;


* limit run time
option ResLim = 100;

* optimality gap
option optCR = .001;

* print solution file
option solprint = on;

* shorten output colums and rows
*option limrow = 0;
*option limcol = 0;

* run in parallel with  threads
option threads = 4;

* use Cplex solver for linear programs
option LP = Cplex;

* use OPT file cplex for minimizing lp problems concurrently
opf.OptFile = 1

option dmpsym;
Solve opf using MIP minimizing z;





Parameters
scen_cost(c,subs,subh)
tot_cost
lines_built(cand,i,j,c)
lb
;

scen_cost(c,subs,subh) = zh.l(c,subs,subh);
tot_cost = z.l;
lines_built(cand,i,j,c)$ii(cand,i,j) = y.l(cand,i,j,c);
lb = opf.objest;
execute_unload 'results', scen_cost, tot_cost, lines_built, lb;
