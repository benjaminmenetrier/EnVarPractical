#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to practice the implementation of EnVar methods.

Goal: implement the different B matrices as full matrices and as sequences of sparse operators, and compare them with Dirac tests

Content:
- Exercise 1: pure sampled covariance matrix
- Exercise 2: localized covariance matrix
- Exercise 3: hybrid covariance matrix

Requires the tools available in "utilities.py".

Programming: Benjamin Menetrier (benjamin.menetrier@met.no)

Licensing: this code is distributed under the APACHE 2 license
Copyright (c) 2023-2024 Meteorologisk Institutt
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

# Maximum ensemble size
Nmax = 1000

# Set random seed (for reproducibility)
np.random.seed(7)

###############################################################################
# Generate true covariance
###############################################################################
Utrue = define_correlation_sqrt(n, Rcor)
Btrue = np.matmul(Utrue, np.transpose(Utrue))
plot_matrix(Btrue, "True covariance")

###############################################################################
# Generate dirac vector
###############################################################################
dirac = np.zeros((n,1))
dirac[int(n/2)] = 1.0
plot_vector(dirac, "Dirac vector")

###############################################################################
# Generate ensemble
###############################################################################
ens = generate_ensemble(n, Utrue, Nmax)
plot_vector(ens[:,0:10], "Ensemble (first 10 members)")


###############################################################################
###############################################################################
# Exercise 1: pure sampled covariance matrix
###############################################################################
###############################################################################
print("Pure sampled covariance matrix:")

# Exercise 2.1: compute ensemble perturbations and compute pure sampled covariance for the whole ensemble

# Remove mean to get ensemble perturbations
mean = np.sum(ens[:,0:Nmax], 1)/Nmax
pert = np.zeros((n, Nmax))
for p in range(Nmax):
    pert[:,p] = ens[:,p]-mean

# Scale ensemble perturbations with 1/sqrt(N-1)
pert *= 1/np.sqrt(Nmax-1)

# Compute pure sampled covariance
Bens = np.zeros((n,n))
for p in range(Nmax):
    Bens += np.matmul(pert[:,p:p+1], np.transpose(pert[:,p:p+1]))
# Could be also done in one line:
# Bens = np.matmul(pert, np.transpose(pert))

# Plot true and pure sampled covariances
plot_matrices(Btrue, "True covariance", Bens, "Pure sampled covariance\n" + str(Nmax) + " members")

# Exercise 1.2: play with the ensemble size N
N_list = [10, 50, 250]
for N in N_list:
    mean = np.sum(ens[:,0:N], 1)/N
    pert = np.zeros((n, N))
    for p in range(N):
        pert[:,p] = ens[:,p]-mean
    pert *= 1/np.sqrt(N-1)
    Bens = np.matmul(pert, np.transpose(pert))
    plot_matrices(Btrue, "True covariance", Bens, "Pure sampled covariance\n" + str(N) + " members")

# Exercise 1.3: compute a pure sampled covariance for 20 members (and don't overwrite it afterwards), apply it to the dirac vector and plot the result
N = 20
mean = np.sum(ens[:,0:N], 1)/N
pert = np.zeros((n, N))
for p in range(N):
    pert[:,p] = ens[:,p]-mean
pert *= 1/np.sqrt(N-1)
Bens = np.matmul(pert, np.transpose(pert))
dirac_test = np.matmul(Bens, dirac)
plot_vector(dirac_test, "Dirac test: pure sampled covariance, 20 members, full matrix")

# Exercise 1.4: dirac test on the pure sampled covariance with 20 members, sequence of operators

# 1.4.1: Apply square-root adjoint
v = np.zeros(N)
for p in range(N):
    v[p] = sum(pert[:,p]*dirac[:,0])

# 1.4.2: Apply square-root
dirac_test_seq = np.zeros((n,1))
for p in range(N):
    dirac_test_seq[:,0] += v[p]*pert[:,p]

# 1.4.3: Plot result
print("-> difference between both formulations: " + str(np.sum(np.abs(dirac_test_seq-dirac_test))))
plot_vector(dirac_test_seq, "Dirac test: pure sampled covariance, 20 members, sequence of operators")


###############################################################################
###############################################################################
# Exercise 2: localized covariance matrix
###############################################################################
###############################################################################
print("Localized covariance matrix:")

# Exercise 2.1: create a localization matrix
Rloc = 25.0
Ul = define_correlation_sqrt(n, Rloc)
L = np.matmul(Ul, np.transpose(Ul))

# Exercise 2.2: apply the localization to the pure sampled covariance (Schur product), plot the localization and the localized covariance
Bloc = Bens*L
plot_matrices(L, "Localization, Rloc = " + str(Rloc), Bloc, "Localized covariance")

# Exercise 2.3: play with the localization radius Rloc
Rloc_list = [50.0, 25.0, 10.0, 1.0]
for Rloc in Rloc_list:
    Ul = define_correlation_sqrt(n, Rloc)
    L = np.matmul(Ul, np.transpose(Ul))
    Bloc = Bens*L
    plot_matrices(L, "Localization, Rloc = " + str(Rloc), Bloc, "Localized covariance")

# Exercise 2.4: dirac test on the localized covariance with Rloc = 25
Rloc = 25.0
Ul = define_correlation_sqrt(n, Rloc)
L = np.matmul(Ul, np.transpose(Ul))
Bloc = Bens*L
dirac_test = np.matmul(Bloc, dirac)
plot_vector(dirac_test, "Dirac test: localized covariance, full matrix")

# Exercise 2.5: dirac test on the localized covariance with Rloc = 25, sequence of operators

# 2.5.1: Apply square-root adjoint
v = np.zeros((n, N))
for p in range(N):
    xtmp = pert[:,p]*dirac[:,0]
    v[:,p] = np.matmul(np.transpose(Ul),xtmp)

# 2.5.2: Apply square-root
dirac_test_seq = np.zeros((n,1))
for p in range(N):
    alpha_p = np.matmul(Ul,v[:,p])
    dirac_test_seq[:,0] += pert[:,p]*alpha_p

 # 2.5.3: Plot result
print("-> difference between both formulations: " + str(np.sum(np.abs(dirac_test_seq-dirac_test))))
plot_vector(dirac_test_seq, "Dirac test: localized covariance, sequence of operators")


###############################################################################
###############################################################################
# Exercise 3: hybrid covariance matrix
###############################################################################
###############################################################################
print("Hybrid covariance matrix:")

# Exercise 3.1: generate localized covariance with localization cut-off length-scale 20.0
Rloc = 20.0
Ul = define_correlation_sqrt(n, Rloc)
L = np.matmul(Ul, np.transpose(Ul))
Bloc = Bens*L

# Exercise 3.2: generate a static B of cut-off length-scale 50.0
Rsta = 50.0
Us = define_correlation_sqrt(n, Rsta)
Bs = np.matmul(Us, np.transpose(Us))
plot_matrix(Bs, "Static covariance")

# Exercise 3.3: define a hybrid covariance matrix with weights 0.7 on the ensemble component and 0.3 on the static component
gammae = 0.7
gammas = 0.3
Bhyb = gammae*Bloc + gammas*Bs
plot_matrix(Bhyb, "Hybrid covariance, weights = " + str(gammae) + " / " + str(gammas))

# Exercise 3.4: play with the hybrid weights, keeping the total variance constant (gammae + gammas = 1)
gammae_list = [0.0, 0.33, 0.66, 1.0]
for gammae in gammae_list:
    gammas = round(1.0-gammae,2)
    Bhyb = gammae*Bloc + gammas*Bs
    plot_matrix(Bhyb, "Hybrid covariance, weights = " + str(gammae) + " / " + str(gammas))

# Exercise 3.5: dirac test on the localized covariance with weights = 0.6 / 0.4, full matrix
gammae = 0.6
gammas = 0.4
Bhyb = gammae*Bloc + gammas*Bs
dirac_test = np.matmul(Bhyb, dirac)
plot_vector(dirac_test, "Dirac test: hybrid covariance, full matrix")

# Exercise 3.6: dirac test on the localized covariance with weights = 0.6 / 0.4, sequence of operators

# 3.6.1: Apply square-root adjoint
v = np.zeros((n, N+1))
betae = np.sqrt(gammae)
betas = np.sqrt(gammas)
for p in range(N):
    xtmp = pert[:,p]*dirac[:,0]
    v[:,p] = betae*np.matmul(np.transpose(Ul),xtmp)
v[:,N] = betas*np.matmul(np.transpose(Us),dirac[:,0])

# 3.6.2: Apply square-root
dirac_test_seq = np.zeros((n,1))
for p in range(N):
    alpha_p = np.matmul(Ul,v[:,p])
    dirac_test_seq[:,0] += betae*pert[:,p]*alpha_p
xtmp = np.matmul(Us,v[:,N])
dirac_test_seq[:,0] += betas*xtmp

# 3.6.3: Plot result
print("-> difference between both formulations: " + str(np.sum(np.abs(dirac_test_seq-dirac_test))))
plot_vector(dirac_test, "Dirac test: localized covariance, sequence of operators")
