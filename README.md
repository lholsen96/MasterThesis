# Master's Thesis by Lars Olsen: "Likelihood-Based Boosting: Approximate Confidence Bands and Intervals for Generalized Additive Models"
R codes used in my master's thesis "Likelihood-Based Boosting: Approximate Confidence Bands and Intervals for Generalized Additive Models" delivered 15th of May 2020 

The main files of interest are:

GAMBoost_stumps: contains GAMBoost with penalized stumps algorithm, alongside all auxiliary functions.

GAMBoost_stumps_with_intercept_update: Same as above, but this implementation of GAMBoost with penalized stumps conduct intercept updates. It relies on auxiliary functions from GAMBoost_stumps.

GAMBoost_splines: contains auxiliary functions to Binder's R-package 'GAMBoost'

GAMBoost_common: contains auxiliary functions that work for both GAMBoost with penalized B-splines and penalized stumps. 

PartBoostR_source_code: contains the implementation of PartBoostR and all auxillary functions.

GenPartBoostR_source_code: contains the implementation of GenPartBoostR and all auxillary functions.

Chapter443_coverage_GAMBoost_splines: contains the code which computes the median coverage of the approximate confidence bands and intervals for f_j and \mu, respectively, in the case of penalized B-splines as base learners. This corresponds to Tables 4.1 - 4.4.

Chapter444_coverage_GAMBoost_stumps: contains the code which computes the median coverage of the approximate confidence bands and intervals for f_j and \mu, respectively, in the case of penalized stumps as base learners. This corresponds to Tables 4.5 - 4.9.

Chapter2_Boostr_PartBoostR_GenPartBoostR_figures: generates the figures in Chapter 2.

Chapter3_Splines_Trees_figures: generates the figures in Chapter 3.

Some of the code sections rely on saved R objects, which was not allowed to upload to GitHub. 
Contact me via my e-mail account to receive these RDS saves. 
