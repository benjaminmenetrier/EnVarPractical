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

# Remove mean to get ensemble perturbations
mean = np.sum(ens[:,0:N], 1)/N
pert = np.zeros((n, N))
for j in range(N):
    pert[:,j] = ens[:,j]-mean

# Scale ensemble perturbations with 1/sqrt(N-1)
pert = pert/np.sqrt(N-1)

# Compute and plot pure sampled covariance
Bens = np.matmul(pert, np.transpose(pert))
plot_matrices(Btar, "Target covariance", Bens, "Pure sampled covariance, " + str(N) + " members")

# Exercise 1.2: play with the ensemble size N
N_list = [10, 50, 250]
for N in N_list:
    mean = np.sum(ens[:,0:N], 1)/N
    pert = np.zeros((n, N))
    for j in range(N):
        pert[:,j] = ens[:,j]-mean
    pert = pert/np.sqrt(N-1)
    Bens = np.matmul(pert, np.transpose(pert))
    plot_matrices(Btar, "Target covariance", Bens, "Pure sampled covariance, " + str(N) + " members")

# Exercise 1.3: dirac test on the pure sampled covariance
dirac_test = np.matmul(Bens, dirac)
plot_vector(dirac_test, "Dirac test: pure sampled covariance, matrix formulation")

# Exercise 1.4: dirac test on the pure sampled covariance, operator formulation

# 1.4.1: Apply square-root adjoint
v = np.zeros(N)
for j in range(N):
    v[j] = sum(pert[:,j]*dirac[:,0])

# 1.4.2: Apply square-root
dirac_test_op = np.zeros((n,1))
for j in range(N):
    dirac_test_op[:,0] = dirac_test_op[:,0] + v[j]*pert[:,j]

# 1.4.3: Plot result
print("Difference between matrix and operator formulations: " + str(np.sum(np.abs(dirac_test_op-dirac_test))))
plot_vector(dirac_test_op, "Dirac test: pure sampled covariance, operator formulation")



###############################################################################
###############################################################################
# Exercise 2: about the localized covariance matrix
###############################################################################
###############################################################################

# Exercise 2.1: create perturbations and pure sampled covariance for 20 members
N = 20
mean = np.sum(ens[:,0:N], 1)/N
pert = np.zeros((n, N))
for j in range(N):
    pert[:,j] = ens[:,j]-mean
pert = pert/np.sqrt(N-1)
Bens = np.matmul(pert, np.transpose(pert))

# Exercise 2.2: create localization and apply it to the pure sampled covariance (Schur product)
Rloc = 25.0
Ul = define_correlation_sqrt(n, Rloc)
L = np.matmul(Ul, np.transpose(Ul))
Bens_L = Bens*L
plot_matrices(L, "Localization, Rloc = " + str(Rloc), Bens_L, "Localized covariance")

# Exercise 2.3: play with the localization radius Rloc
Rloc_list = [50.0, 25.0, 10.0, 1.0]
for Rloc in Rloc_list:
    Ul = define_correlation_sqrt(n, Rloc)
    L = np.matmul(Ul, np.transpose(Ul))
    Bens_L = Bens*L
    plot_matrices(L, "Localization, Rloc = " + str(Rloc), Bens_L, "Localized covariance")

# Exercise 2.4: dirac test on the localized covariance with Rloc = 25
Rloc = 25.0
Ul = define_correlation_sqrt(n, Rloc)
L = np.matmul(Ul, np.transpose(Ul))
Bens_L = Bens*L
dirac_test = np.matmul(Bens_L, dirac)
plot_vector(dirac_test, "Dirac test: localized covariance, matrix formulation")

# Exercise 2.5: dirac test on the localized covariance with Rloc = 25, operator formulation

# 2.5.1: Apply square-root adjoint
v = np.zeros((n, N))
for j in range(N):
    xtmp = pert[:,j]*dirac[:,0]
    v[:,j] = np.matmul(np.transpose(Ul),xtmp)

# 2.5.2: Apply square-root
dirac_test = np.zeros((n,1))
for j in range(N):
    xtmp = np.matmul(Ul,v[:,j])
    dirac_test[:,0] = dirac_test[:,0] + pert[:,j]*xtmp

 # 2.5.3:Plot result
plot_vector(dirac_test, "Dirac test: localized covariance, operator formulation")



###############################################################################
###############################################################################
# Exercise 3: about the hybrid covariance matrix
###############################################################################
###############################################################################

# Exercise 3.1: generate localized covariance with localization cut-off length-scale 20.0
Rloc = 20.0
Ul = define_correlation_sqrt(n, Rloc)
L = np.matmul(Ul, np.transpose(Ul))
Bens_L = Bens*L

# Exercise 3.2: generate a static B of cut-off length-scale 50.0
Rsta = 50.0
Usta = define_correlation_sqrt(n, Rsta)
Bsta = np.matmul(Usta, np.transpose(Usta))
plot_matrix(Bsta, "Static covariance")

# Exercise 3.3: define a hybrid covariance matrix with weights 0.7 on the ensemble component and 0.3 on the static component
beta_ens = 0.7
beta_sta = 0.3
Bhyb = beta_ens*Bens_L + beta_sta*Bsta
plot_matrix(Bhyb, "Hybrid covariance, weights = " + str(beta_ens) + " / " + str(beta_sta))

# Exercise 3.4: play with the hybrid weights, keeping the total variance constant
beta_ens_list = [0.0, 0.33, 0.66, 1.0]
for beta_ens in beta_ens_list:
    beta_sta = round(1.0-beta_ens,2)
    Bhyb = beta_ens*Bens_L + beta_sta*Bsta
    plot_matrix(Bhyb, "Hybrid covariance, weights = " + str(beta_ens) + " / " + str(beta_sta))

# Exercise 3.5: dirac test on the localized covariance with weights = 0.6 / 0.4, matrix formulation
beta_ens = 0.6
beta_sta = 0.4
Bhyb = beta_ens*Bens_L + beta_sta*Bsta
dirac_test = np.matmul(Bhyb, dirac)
plot_vector(dirac_test, "Dirac test: hybrid covariance, matrix formulation")

# Exercise 3.6: dirac test on the localized covariance with weights = 0.6 / 0.4, operator formulation

# 3.6.1: Apply square-root adjoint
v = np.zeros((n, N+1))
for j in range(N):
    xtmp = np.sqrt(beta_ens)*dirac[:,0]
    xtmp = pert[:,j]*xtmp
    v[:,j] = np.matmul(np.transpose(Ul),xtmp)
xtmp = np.sqrt(beta_sta)*dirac[:,0]
v[:,N] = np.matmul(np.transpose(Usta),xtmp)

# 3.6.2: Apply square-root
dirac_test_op = np.zeros((n,1))
for j in range(N):
    xtmp = np.matmul(Ul,v[:,j])
    xtmp = pert[:,j]*xtmp
    dirac_test_op[:,0] = dirac_test_op[:,0] + np.sqrt(beta_ens)*xtmp
xtmp = np.matmul(Usta,v[:,N])
dirac_test_op[:,0] = dirac_test_op[:,0] + np.sqrt(beta_sta)*xtmp

# 3.6.3: Plot result
print("Difference between matrix and operator formulations: " + str(np.sum(np.abs(dirac_test_op-dirac_test))))
plot_vector(dirac_test, "Dirac test: localized covariance, operator formulation")
