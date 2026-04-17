# Nonlinear Bayesian filter with uncertain loads

This repository implements the Nonlinear Bayesian filter for structural dynamics with uncertain forces [CITE PENDING]. This filter estimates states and uncertain forces using Gaussian random walks, the unscented transform, and single-step time-integration algorithms.

The repository contains a [numerical validation and comparison with the UKF](./numerical_validation.m), a [sensitivity analysis on the noise level](./sensitivity_measurement_noise.m), and a [sensitivity analysis on the number of measurements](./sensitivity_number_measurement.m). In the examples, the analysed system is a 6-DOFs shear chain with hysteretic springs governed by the Bouc-Wen model.
