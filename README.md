# Nonlinear Bayesian filter with uncertain forces

This repository implements the Nonlinear Bayesian filter for structural dynamics with uncertain forces [CITE PENDING]. This filter estimates states and uncertain forces using Gaussian random walks, the unscented transfrom, and single-step time-integration algorithms.

The repository contains a numerical validation and comparison with the UKF, a sensitivity analysis for the noise level, and a sensitivity analysis with the number of measurements. In the examples, the analysed system is a 6-DOFs shear chain with hysteretic springs governed by the Bouc-Wen model.
