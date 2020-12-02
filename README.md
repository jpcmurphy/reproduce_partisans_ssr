# partisans_revisit

Code to reproduce analyses in "[Issue Alignment and Partisanship in the American Public: Revisiting the 'Partisans without Constraint' Thesis](https://osf.io/preprints/socarxiv/jex9k/)." Scripts should be run in following order.

(1) __makeBaldGelData.R__      Makes main data set of inter-item correlations by year including item-domain indicators <br/>
(2) __mixedFXsmoothsplines.R__      Mixed effects smoothing spline model for aggregate trend and graphs results <br/>
(3) __gamm_by_domain.R__      Generalized additive mixed effects models for trends by issue domain and graphs results <br/>
(4) __subgroupData.R__      Calculates inter-item correlations by demographic subgroup <br/>
(5) __subgroupGAMM.R__      Runs GAMMs by issue domain for each subgroup and graphs results <br/>
(6) __media_boxes16.R__     Box plots of polarization measures by media engagement and consumption <br/>

<br/>

Versions of libaries used in model fitting: sme (v1.0.2); gamm4 (v0.2-6); mgcv (v1.8-28); lme4 (v1.1-23). <br/>

The scripts reference a top-level (i.e. this level) folder _data/_, which should store both the raw ANES election file input and the SQLite database created/modified in __makeBaldGelData.R__ and __subgroupData.R__. <br/>

The analyses rely on the [American National Election Studies Cumulative Data File (1948-2016)](https://electionstudies.org/data-center/anes-time-series-cumulative-data-file/) as of December 6, 2018. (ANES_CDF_VERSION:2018-Dec-06)
