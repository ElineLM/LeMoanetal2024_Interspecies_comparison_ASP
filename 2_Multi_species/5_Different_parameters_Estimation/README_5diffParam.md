These files try to estimate parameter values for the 5 species at the same time, with some parameters being equal, other being forced by physical co-variation rules and other being different but with a weight to control how much they can be different between species.

In this study, the last possibility was tried for pM and v parameters.

For each estimation, it is needed to
1. modify the pars_init_group, by modifying the weight of one or two parameters
2. modify initial values of the parameters, if we want values to be different between species, we need to give 5 initial values (5 being the number of species estimated)
3. In the run_group file, do not forget to modify the name of the file in which results will be saved to not erase the previous simulations

In the study, estimations were realised with 3 runs, as described in the run file
