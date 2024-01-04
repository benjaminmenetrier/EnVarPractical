#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to practice the implementation of EnVar methods:
- Exercise 1: about the pure sampled covariance matrix
- Exercise 2: about the localized covariance matrix
- Exercise 3: about the hybrid covariance matrix

Requires the tools available in "utilities.py".

Data assimilation formation at Cerfacs 2023.

Programming: Benjamin Menetrier, 2023 (benjamin.menetrier@met.no)

Licensing: this code is distributed under the APACHE 2 license
  Copyright (c) 2023 Meteorologisk Institutt
"""
###############################################################################
# Imports
###############################################################################
import numpy as np
from utilities import *

###############################################################################
# General parameters (should not be changed)
###############################################################################

# Domain size
n = 101

# Correlation cut-off length-scale (unit: number of points)
Rcor = 15.0

# Ensemble size
N = 1000

# Set random seed (for reproducibility)
np.random.seed(7)

###############################################################################
# Generate target covariance
###############################################################################
Utar = define_correlation_sqrt(n, Rcor)
Btar = np.matmul(Utar, np.transpose(Utar))
plot_matrix(Btar, "Target covariance")

###############################################################################
# Generate dirac vector
###############################################################################
dirac = np.zeros((n,1))
dirac[int(n/2)] = 1.0
plot_vector(dirac, "Dirac vector")

###############################################################################
# Generate ensemble
###############################################################################
ens = generate_ensemble(n, Rcor, N)
plot_vector(ens[:,0:10], "Ensemble (first 10 members)")



###############################################################################
###############################################################################
# Exercise 1: about the pure sampled covariance matrix
###############################################################################
###############################################################################

# Exercise 1.1: compute the pure sampled covariance using ensemble "ens" (shape n x N)
# ...

# Exercise 1.2: play with the ensemble size N
N_list = [10, 50, 250]
for N in N_list:
    # ...

# Exercise 1.3: dirac test on the pure sampled covariance
# ...

# Exercise 1.4: dirac test on the pure sampled covariance, operator formulation

# Apply square-root adjoint
# ...

# Apply square-root
# ...

# Plot result
# ...



###############################################################################
###############################################################################
# Exercise 2: about the localized covariance matrix
###############################################################################
###############################################################################

# Exercise 2.1: create perturbations and pure sampled covariance for 20 members
# ...

# Exercise 2.2: create localization and apply it to the pure sampled covariance (Schur product)
# ...

# Exercise 2.3: play with the localization radius Rloc
Rloc_list = [50.0, 25.0, 10.0, 1.0]
for Rloc in Rloc_list:
    # ...

# Exercise 2.4: dirac test on the localized covariance with Rloc = 25
# ...

# Exercise 2.5: dirac test on the localized covariance with Rloc = 25, operator formulation

# Apply square-root adjoint
# ...

# Apply square-root
# ...

# Plot result
# ...



###############################################################################
###############################################################################
# Exercise 3: about the hybrid covariance matrix
###############################################################################
###############################################################################

# Exercise 3.1: generate localized covariance with localization cut-off length-scale 20.0
# ...

# Exercise 3.2: generate a static B of cut-off length-scale 50.0
# ...

# Exercise 3.3: define a hybrid covariance matrix with weights 0.7 on the ensemble component and 0.3 on the static component
# ...

# Exercise 3.4: play with the hybrid weights, keeping the total variance constant
beta_ens_list = [0.0, 0.33, 0.66, 1.0]
for beta_ens in beta_ens_list:
    # ...

# Exercise 3.5: dirac test on the localized covariance with weights = 0.6 / 0.4, matrix formulation
# ...

# Exercise 3.6: dirac test on the localized covariance with weights = 0.6 / 0.4, operator formulation

# Apply square-root adjoint
# ...

# Apply square-root
# ...

# Plot result
# ...
