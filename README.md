OpenFOAM-AD
===========

[![Regression Test](https://github.com/DAFoam/OpenFOAM-AD/actions/workflows/reg_tests.yml/badge.svg)](https://github.com/DAFoam/OpenFOAM-AD/actions/workflows/reg_tests.yml)

This repository contains the OpenFOAM source codes differentiated by automatic differentiation (AD) in forward and reverse modes. OpenFOAM is a free, open source CFD software [released and developed primarily by OpenCFD Ltd](http://www.openfoam.com).

Download
--------

We only differentiate certain versions of OpenFOAM releases. So check the branches (e.g., openfoam-v1812-ad) on this repo; **the main branch does not contain any AD codes**.

Installation
------------

The installation of OpenFOAM-AD is similar to that of OpenFOAM. One needs to first install all prerequisites, source the etc/bashrc file, and then run `./Allwmake`.

The default build will be for forward mode AD. To compile reverse mode AD, change `WM_CODI_AD_MODE` to `CODI_AD_REVERSE` in etc/bashrc, source it, and rebuild.

NOTE: OpenFOAM-AD only differentiates necessary libraries for computing partial derivatives and matrix-vector products for [DAFoam](https://dafoam.github.io), it has NOT differentiated the entire OpenFOAM code yet. In other words, some functionalities are still missing (e.g. combustion models).

Acknowledgement
---------------

OpenFOAM-AD uses the AD tools [CoDiPack](https://github.com/scicompkl/codipack) and [MeDiPack](https://github.com/scicompkl/medipack), developed by Dr. Nicolas Gauger's group at TU Kaiserslautern. In addition, the differentiation is inspired by [discreteAdjointOpenFOAM](https://www.stce.rwth-aachen.de/research/software/discreteadjointopenfoam), developed by Dr. Markus Towara et al. from RWTH Aachen University.

License
-------

Copyright 2021 iDesign Lab, Aerospace Engineering Department, Iowa State University.

Distributed using the GNU General Public License (GPL), version 3; see the COPYING file for details.
