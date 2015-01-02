Test Problems on a Rectangular Domain
=====================================

The scripts in this directory solve three model problems on
a rectangular domain.  In each case, a simple analytical solution
is known, and used to compute the maximum nodal error in the finite
element solution.  A sequence of mesh refinements allows us to estimate
the convergence rate k such that the error is of order h^k.  We expect
k=2.

First run the script setup.jl which generates the sequence of mesh files.
Then run individually

dirichlet.jl
neumann.jl
poisson.jl

See the leading comments in each file for a precise mathematical
description of the boundary value problem.
