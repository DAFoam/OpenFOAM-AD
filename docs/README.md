OpenFOAM-AD
===========

This repository contains the OpenFOAM source code, differentiated using automatic differentiation (AD) in both forward and reverse modes. OpenFOAM is a free, open source CFD software [released and developed primarily by OpenCFD Ltd](http://www.openfoam.com).

Download
--------

We only differentiate certain versions of OpenFOAM releases. Please check the branches (e.g., v2506-ad) in this repository; **the main branch does not contain any AD code**.

Installation
------------

The default build will be for the reverse mode AD (`export WM_AD_MODE=ADR`). To compile the forward mode AD, change `WM_AD_MODE` to `ADF` in OpenFOAM-AD/etc/bashrc, source it, and rebuild.

NOTE: OpenFOAM-AD only differentiates necessary libraries for computing partial derivatives and matrix-vector products for [DAFoam](https://dafoam.github.io); it has NOT differentiated the entire OpenFOAM code yet. In other words, some functionalities are still missing (e.g., functionObjects and Thirdparty-related libs).

Acknowledgement
---------------

OpenFOAM-AD uses the AD tools [CoDiPack](https://github.com/scicompkl/codipack) and [MeDiPack](https://github.com/scicompkl/medipack), developed by Dr. Nicolas Gauger's group at TU Kaiserslautern. In addition, the differentiation is inspired by [discreteAdjointOpenFOAM](https://gitlab.stce.rwth-aachen.de/towara/discreteadjointopenfoam_adwrapper), developed by Dr. Markus Towara et al. from RWTH Aachen University.

License
-------

Distributed using the GNU General Public License (GPL), version 3; see the COPYING file for details.
