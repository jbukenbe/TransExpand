function problem = run_tep2_gams(problem)
% this function takes the problem input for the transmission expansion
% problem and formulates it into a GAMS MIP that solves the two stage
% stochastic TEP. Clustering information for the second stage decision
% nodes is needed for compressed models.


%History            
%Version    Date        Who     Summary
%1          05/18/2018  JesseB  Adapted from run_tep_gams