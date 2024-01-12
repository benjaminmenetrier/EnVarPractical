#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tools required for the EnVar exercises:
- define_correlation_sqrt(n, R): return the square-root of a correlation matrix of size (n,n), with a correlation cut-off length-scale R
- generate_ensemble: returns an ensemble of size nens, randomizing the correlation matrix of size (n,n), with a correlation cut-off length-scale R
- plot_vector(X, title): plot a vector (or a set of vectors) X with the provided title
- plot_matrix(A, title): plot matrix A with the provided title
- plot_matrices(A1, title1, A2, title2): plot matrix A1 and A2 with the provided titles title1 and title2, respectively.

Programming: Benjamin Menetrier (benjamin.menetrier@met.no)

Licensing: this code is distributed under the APACHE 2 license
Copyright (c) 2023-2024 Meteorologisk Institutt
"""
###############################################################################
# Imports
###############################################################################
import matplotlib.pyplot as plt
import numpy as np

###############################################################################
# define_correlation_sqrt
###############################################################################
def define_correlation_sqrt(n, R):
    U = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            dist = np.abs(i-j)
            if (dist>n/2):
                dist = n-dist
            dist = dist/R
            if dist <= 0.5:
                U[i,j] = 1.0-(2.0*dist)
    x = np.zeros((n, 1))
    x[0,0] = 1.0
    x = np.matmul(np.transpose(U), x)
    x = np.matmul(U, x)
    U = U/np.sqrt(x[0,0])
    return U

###############################################################################
# generate_ensemble
###############################################################################
def generate_ensemble(n, Ul, nens):
    # Generate Gaussian deviates
    X = np.random.normal(0.0, 1.0, (n, nens))

    # Apply localization square-root
    XX = np.zeros((n, nens))
    for j in range(nens):
        XX[:,j] = np.matmul(Ul, X[:,j])

    return XX

###############################################################################
# plot_vector
###############################################################################
def plot_vector(X, title):
    n, nens = np.shape(X)
    fig, ax = plt.subplots(figsize=(14,7))
    ax.set_xlim([0,n-1])
    ax.set_title(title)
    ax.axhline(y=0, color='r')
    ax.plot(X, 'k')
    plt.show()

###############################################################################
# plot_matrix
###############################################################################
def plot_matrix(A, title):
    fig = plt.figure()
    plt.set_cmap('jet')
    ax = fig.add_subplot(111)
    ax.set_title(title)
    cax = ax.matshow(A, vmin=-0.2, vmax=1.2)
    fig.colorbar(cax)
    plt.show()

###############################################################################
# plot_matrices
###############################################################################
def plot_matrices(A1, title1, A2, title2):
    fig = plt.figure()
    plt.set_cmap('jet')
    ax = fig.add_subplot(121)
    ax.set_title(title1)
    cax = ax.matshow(A1, vmin=-0.2, vmax=1.2)
    ax = fig.add_subplot(122)
    ax.set_title(title2)
    cax = ax.matshow(A2, vmin=-0.2, vmax=1.2)
    plt.show()

