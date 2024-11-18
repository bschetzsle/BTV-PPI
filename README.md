# BTV-PPI

Bayesian Time-varying Psychophysiological Interaction Model

This repository contains the code for applying our proposed Bayesian Time-varying Psychophysiological Interaction Model to both simulated data and real-world data.

The directory "Simulation Code" contains R scripts to simulate data and then fit both the standard Generalized PPI model and our proposed extension to this simulated data. Partial correlations between simulated brain regions are derived from the fitted model using "calculate_par_cors.R" script. The mean squared error of these estimated partial correlations can be calculated using "get_MSE_from_par_cors.R". There are additional files that generate plots used in the paper.

The directory "Applied Code" contains the code used to fit our proposed model to real-world data, though does not contain that data. There are several files that fit our model to the data which consider different formulations of the stimulus. The R script "fit_entropy_difference_separate_stimulus_effects.R" considers two stimuli: a continuous measure of difference in entropy between the fast and slow learning rates, and an indicator for self-self transitions, when the subject was shown the same image twice in a row. "fit_entropy_difference.R" does not include self-self transitions as a stimulus. "fit_entropy_separate_stimulus_effects.R" considers four different stimuli: image presentation, button press, entropy under a slow learning rate and entropy under a fast learning rate. These scripts all return the matrix describing the estimated partial correlations between pairs of brain regions. The script "perform_selection_on_btv_ppi_partial_correlations.R" and "perform_selection_on_btv_ppi_partial_correlations_entropy_difference.R" perform the selection method of Chandra and Bhattacharya (2019) on the estimated partial correlations. There are some additional files the visualize the results for use in the paper.
