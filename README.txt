This repository stores experiments for accelerating the solve time for the transmission expansion problem with multivariate information obtained through the singular value decomposition. 

** For Mort **
The code can be tested by running: 
>> dyn_pls_demo()
Parts of the code can be run in parallel so running >>parpool(‘local’) fist will significantly speed up execution.



Changes to the model parameters can be made in the DC_OPF_init.m file. The parameters that will make the most impact are:
params.svd.scen_n 
This changes how many scenarios are used in the latent factor approximation after the initial samples are finished running.

params.initial_samp_n 
This changes the number of plans that are sampled and used for the SVD (only in the first step of the algorithm) **NOTE: the function >> fraction_fact_samp.m in the Support Functions folder may ignore this value if it is too small

params.refine_samp_n
This changes the number of plans that are sampled in all steps after the SVD has been performed (steps 2 – Z) **NOTE: the function >> fraction_fact_samp.m in the Support Functions folder may ignore this value if it is too small

params.maxz
This changes the number of iterations the algorithm will perform before terminating.

A larger set of scenarios can be used by changing the input data in the %% Read GAMS Data section to include the commented out sections.




To experiment with using the results from the SVD matrix, the parts of the code that are most relevant are as follows:

1) In dyn_pls_demo.m, the part commented with %% Save SVD Latent Factor Info
This part of the code performs SVD (the U matrix is not currently saved) and selects the most significant scenarios to examine by selecting the scenarios with large values in the V matrix. In the future this will be replaced with a determinant maximizing function based on Fedorov’s exchange algorithm. 
The number of scenarios chosen is params.svd.scen_n. The important scenarios are stored as latent_scen while the other scenarios are stored as filler_scen. The best candidate solution value is also saved and if the approximate value of a new candidate plan is better than the best known plan, the filler_scen will be calculated as well for precision.

2) In par_plan_ops.m the part commented with % estimate operations with svd latent factor approximation
This section runs the scenarios marked with latent_scen and then calculates the approximate cost of all the other scenarios, if the approximate cost of the plan is better than the designated threshold, the filler_scen are calculated to confirm the true cost of the plan. If the approximate cost is worse than the best known plan, the other scenarios are skipped.

3) mat_fill_svd.m in support functions
This performs the latent factor approximation of an incomplete input vector. This is done with gradient descent here, but the normal equations can be solved directly so this will be replaced with the material in item 4 soon.

4) svd_approx.m in Latent Factor Test Functions
This function isolates most of the SVD approximation components for testing. The main result is that the hat matrix H can be computed directly and multiplied by the partial input vector to produce the approximate output instead of finding it with gradient descent.




Ideas for future work
Information about each plan can be gained from saving the U matrix from the SVD. A new candidate plan can be found with a set cover type problem of the form:

Variable Def:
Uxk = estimated U vector for new composite plan (made by combining known plans)
pi = proportion of known plan i to include in the new composite plan

Objective function:
Min z = Ux*S*V'

s.t.

Uxk = sum(i, pi*Uik) for all k		Approximate U vector is combination of U vecors of plans included in composite plan

-minj(Ui.k)*c < Ux < maxj(Ui.k)*c	Don't let Ux get too far away from what is known to be possible for each latent factor 

sum(i, pi*cand_plan_has_line_l) < 1 for all l 	Can't build same line multiple times through different plans, each line can only be built once

pi > 0 for all i			No negative proportions 
