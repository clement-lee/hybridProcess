<!-- README.md is generated from README.Rmd. Please edit that file -->
hybridProcess
=============

The goal of hybridProcess is to provide functions for fitting and
simulation of hybrid processes, one kind of non-homogeneous Poisson
processes (NHPP). Specifically, the intensity of this kind of NHPP is
proportional to the product of a power term, and an exponential term, of
time.

Installation
------------

You can install hybridProcess from github with:

    # install.packages("devtools")
    devtools::install_github("clement-lee/hybridProcess")

Example
-------

Here is an example of simulating from, and fitting, a hybrid process,
with power parameter *λ* = 0.66, exponent *θ* = 0.1, and scale parameter
*ϕ* = 5.0.

    # library(hybridProcess)
    ## Simulation
    set.seed(1234L)
    x <- sim_hybrid(
      T = 1000, x = 0.0, beta = 0.0, kappa = 0.0,
      lambda = 0.66, phi = 5.0, theta = 0.1,
      psi = 0.0, tau = Inf
    )
    ## Fitting
    obj0 <- optim(c(0.1, 0.1, 0.1), llik_nhpp_hybrid, x = x, control = optimctrl())
    obj0$par # compare with true parameter values
    #> [1] 0.6945695 0.1081603 4.0358572

Further information
-------------------

The example above uses a function `sim_hybrid()` that is more flexible
than a simple hybrid process. A wrapper function `sim_wrapper()` enables
simulation of a collection of hybrid processes, each of which can be
initiated at a different time point, with the possible inclusion of a
covariate. Such data can be fitted by a hierarchical model of hybrid
processes, via the function `run_mcmc()`.

For details of the hierarchical model, as well as the definitions of
other parameters in `sim_hybrid()`, please refer to [this ArXiv
paper](https://arxiv.org/abs/1802.01987).
