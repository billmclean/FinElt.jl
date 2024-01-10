FinElt.jl
=========

**This package is no longer maintained.  You might be interested in
[this newer replacement](https://github.com/billmclean/SimpleFiniteElements.jl).**

A simple package for a first course in finite element methods.

These Julia routines rely on [Gmsh][1] for mesh generation and
visualization.  Please see the Quick Start Guide (under doc/) for
examples illustrating the use of FinElt.

Note that FinElt supports only msh file format version 2.2 so you need
to use the `-format msh22` command-line option for `gmsh`.

Changelog:

13-01-2020: v0.5.2 Fixes to ensure use of Gmsh 2.2 msh file format.

20-02-2017: v0.0.3 Add support for variable coefficient.

16-07-2016: v0.0.2 Change ASCIIString to String for Julia 0.5.

26-09-2015: v0.0.1

[1]: http://geuz.org/gmsh/

[![Build Status](https://travis-ci.org/billmclean/FinElt.jl.svg?branch=master)](https://travis-ci.org/billmclean/FinElt.jl)
