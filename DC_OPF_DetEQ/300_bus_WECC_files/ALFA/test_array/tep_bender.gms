$ontext
This model solves the transmission expansion plan that minimizes the cost of the
DC OPF model for the WECC 300-Bus test case.

Options:
         Active Scenarios: Choose which scenarios to include in analysis e.g. /s2/ or /s1*s3/ etc...


History
Version__Date____________Who_____Summary
1        06/11/2018      JesseB  copy of lf_tep_bender to run simple weighted avg monte carlo or k means opt




$offtext
$offlisting
$include        "WECC_bender_gams.txt"

SCALARS
ThetaMax         Maximum value of theta [radians]   /.5/
cpns             Cost of non-served power           /1000/
csg              Cost of shed renewable gen         /300/
lb               Bender's lower bound       /-INF/
ub               Bender's upper bound       /INF/
done             Bender's termination indicator     /0/
gap              Bender's optimality gap
;

VARIABLES
zMst                  Value of master problem objective function
z                     Value of the opf objective function
zs(subs)              value of the objective function in each scenario
p(l, i, j, subs)      Active power Flow in line i - j [MW]
theta(i, subs)        Angle of voltage in node i [rad]
;

POSITIVE VARIABLES
rec                    Estimation of the recourse funcion in the master problem
x(g, subs)             Active power dispatched by generator g [MW]
pns(i, subs)           Non-served power in node i [MW]
sg(g, subs)            Generation shed at generator g [MW]
*loss(l,i,j, subs)      Losses in the line i-j [MW]
;

BINARY VARIABLES
y(cand,i,j)            1 if line is built from i to j 0 otherwise
;

EQUATIONS
eMstFun                     Master problem objective function
eCutFun(iter)               Master problem cut

eObFun                          Subproblem objective Function
eSubFun(subs)                   Objective function at each scenario
eDemand(i, subs)                    Demand met from generation and power inflow
ePowerFlow(built,i,j, subs)             DC power flow in existing lines
eCandLineFlowHi(cand,i,j,subs)         DC power flow upper bound in candidate lines
eCandLineFlowLow(cand,i,j,subs)        DC power flow lower bound in candidate lines
eCandCapHi(cand,i,j,subs)              Line capacity upper bound in candidate lines
eCandCapLow(cand,i,j,subs)             Line capacity lower bound in candidate lines
;
eMstFun..                     zMst =e= sum((cand,i,j)$ii(cand,i,j), dtl(cand,i,j,'cost')*y(cand,i,j)) + rec;

eCutFun(dyniter)..            rec =g= pBender_constant(dyniter) + sum((cand,i,j)$ii(cand,i,j), pBender_slope(dyniter,cand,i,j) * y(cand,i,j));


eObFun..                               z =e= sum(subs, zs(subs));

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
eCandLineFlowHi(cand,i,j,subs)$ii(cand,i,j)..         -p(cand,i,j, subs) + 1000*dtl(cand,i,j,'s')*(theta(i, subs) - theta(j, subs)) =g= 10000*(-1+y.l(cand,i,j));
eCandLineFlowLow(cand,i,j,subs)$ii(cand,i,j)..        p(cand,i,j, subs) - 1000*dtl(cand,i,j,'s')*(theta(i, subs) - theta(j, subs)) =g= 10000*(-1+y.l(cand,i,j));

eCandCapHi(cand,i,j,subs)$ii(cand,i,j)..              p(cand,i,j,subs) =l= y.l(cand,i,j)*dtl(cand,i,j,'pmax');
eCandCapLow(cand,i,j,subs)$ii(cand,i,j)..             p(cand,i,j,subs) =g= -y.l(cand,i,j)*dtl(cand,i,j,'pmax');




ii(l,i,j) $ dtl(l,i,j,'pmax') = YES;
x.up(g, subs) = gdata(g,'pmax');
x.lo(g, subs) = gdata(g,'pmin');
p.up(built,i,j, subs)$ii(built,i,j) = dtl(built,i,j,'pmax');
p.lo(built,i,j, subs)$ii(built,i,j) = -dtl(built,i,j,'pmax');
theta.up(i, subs) = thetamax;
theta.lo(i, subs) = -thetamax;

x.fx(nd, subs) = pRenew(nd,subs)*gdata(nd,'pmax');
theta.fx(i, subs) $ [ORD(i) = 1] = 0;


Model master_prob /eMSTFun,eCutFun/
Model opf /eObFun, eSubFun, eDemand, ePowerFlow, eCandCapHi, eCandCapLow, eCandLineFlowHi, eCandLineFlowLow/;


* limit run time
option ResLim = 1200;

* optimality gap
opf.optCR = .0001;
master_prob.optCR = .0001

* dont print solution file
option solprint = off;

* shorten output colums and rows
option limrow = 0;
option limcol = 0;

* run in parallel with threads
option threads = 10;

* use Cplex solver for linear programs
option LP = Cplex;

* use OPT file cplex for minimizing lp problems concurrently
opf.OptFile = 1;

*option dmpsym;




*-------------------- Benders Algorithm --------------------

* Initial master problem plan without cuts
dyniter(iter) = NO;
pBender_constant(iter) = 0;
pBender_slope(iter,cand, i, j) = 0;
Solve master_prob using MIP minimizing zMst;
lb = zMst.l


* Iterate through sub-problem master problem refinement
Loop (iter$(not done),

      dyniter(iter) = yes;

* Solve sub-problems
      solve opf minimizing z using lp;

* Generate new cut with weights
      pBender_constant(iter) = sum(subs, psub_w(subs)*(zs.l(subs)
                       - sum((cand,i,j)$ii(cand,i,j),10000*y.l(cand,i,j)*(eCandLineFlowHi.m(cand,i,j,subs) +  eCandLineFlowLow.m(cand,i,j,subs)))
                       - sum((cand,i,j)$ii(cand,i,j), dtl(cand,i,j,'pmax')*y.l(cand,i,j)*(eCandCapHi.m(cand,i,j,subs)- eCandCapLow.m(cand,i,j,subs)))))/8760;

* abs(p.l(cand,i,j,subs))

      pBender_slope(iter,cand, i, j)$(ii(cand,i,j)) = sum(subs,
                       psub_w(subs)*(10000*(eCandLineFlowHi.m(cand,i,j,subs) + eCandLineFlowLow.m(cand,i,j,subs))
                       + dtl(cand,i,j,'pmax')*(eCandCapHi.m(cand,i,j,subs)- eCandCapLow.m(cand,i,j,subs))))/8760;



*      pBender_constant(iter) =  (sum(subs, psub_w(subs)*(zs.l(subs)
*              - sum((cand,i,j)$ii(cand,i,j),((eCandCapHi.m(cand,i,j,subs)
*              - eCandCapLow.m(cand,i,j,subs))*dtl(cand,i,j,'pmax')*y.l(cand,i,j))+
*              (eCandLineFlowHi.m(cand,i,j,subs) - eCandLineFlowLow.m(cand,i,j,subs))*10000*(1-y.l(cand,i,j)))))/8760;

*      pBender_slope(iter,cand, i, j)$(ii(cand,i,j)) = (sum(subs, psub_w(subs)*((eCandCapHi.m(cand,i,j,subs)
*              - eCandCapLow.m(cand,i,j,subs))*dtl(cand,i,j,'pmax')+
*              (-eCandLineFlowHi.m(cand,i,j,subs)+ eCandLineFlowLow.m(cand,i,j,subs))*10000)))/8760;


* Update upper bound
* ub = min(ub, zMst.l - rec.l + (sum(subs,psub_w(subs)*(zs.l(subs)-psub_mean(subs)))+53000000000)/8760);
      if(ub > (zMst.l - rec.l + sum(subs, psub_w(subs)*zs.l(subs))/8760),
         ub = (zMst.l - rec.l + sum(subs, psub_w(subs)*zs.l(subs))/8760);
         py_best(cand,i,j)$ii(cand,i,j) = y.l(cand,i,j);
         display y.l;
      );
      gap = (ub - lb) / (1 + abs(lb));
      Display ub,lb,gap;

* Test for convergence
      if (gap<0.0001,
          Display "converged";
          if (gap<-0.0001,
              display "with error";
          );
          done = 1;
      else
* Continue iterating, solve master problem with new cuts

      Solve master_prob using MIP minimizing zMst;
      display rec.l;
      lb = zMst.l;
      );
);



Parameters
scen_cost(subs)
tot_cost
lines_built(cand,i,j)
;

scen_cost(subs) = zs.l(subs);
tot_cost = ub;
lines_built(cand,i,j)$ii(cand,i,j) = py_best(cand,i,j);
execute_unload 'results', scen_cost, lines_built, tot_cost, lb;
