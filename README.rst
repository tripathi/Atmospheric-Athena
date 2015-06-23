Atmospheric Athena:
==================
3D Atmospheric escape model with ionizing radiative transfer

What is Atmospheric Athena?
  Atmospheric Athena is a code intended to simulate hydrodynamic escape from close-in giant planets in 3D.  It uses the `Athena hydrodynamics code <https://trac.princeton.edu/Athena/>`_ (v4.1) with a new ionizing radiative transfer implementation based on `Krumholz et al, 2007 <http://arxiv.org/abs/astro-ph/0606539>`_, to self-consistently model photoionization driven winds from the planet.  The code is fully compatible with static mesh refinement and MPI parallelization.  It can handle arbitrary planet potentials and stellar initial conditions.  The physical system that the code models and the equations that it solves are fully described in Tripathi, Kratter, Murray-Clay, & Krumholz, 2015.

Files of note in this repository
  * Mass loss input file: *tst/massloss/athinput.ioniz_sphere_hires*
  * Mass loss problem file: *src/prob/ioniz_sphere.c*
  * Ionization front test input file: *tst/ionradiation/athinput.ifront*
  * Ionization front problem file: *src/prob/ifront.c*

Detailed Doxygen documentation can be found in doc/doxygen/html/index.html.

Configuring and Running
-----------------------
Adapted from the `Athena Tutorial <https://trac.princeton.edu/Athena/wiki/AthenaDocsTut>`_

1. Clean up files from the last compilation.
::
  make clean
2. Configure the file with desired options.  For our 3D mass loss simulations, we used:
::
  ./configure --with-problem=ioniz_sphere --with-gas=hydro --enable-ion-radiation --enable-ion-plane --with-flux=roe --enable-mpi --enable-h-correction --enable-smr
3. Compile. This will create the /bin directory, if it does not already exist. 
::
  make all
4. Run, using the appropriate input file.  If using MPI, use ``mpirun`` and specify the number of processors
::
  bin/athena -i tst/massloss/athinput.ioniz_sphere_hires

For more detailed running instructions, refer to the `Athena documentation <https://trac.princeton.edu/Athena/wiki/AthenaDocs>`_.
  
Copyright
---------
This work is a modified version of `Athena <https://trac.princeton.edu/Athena/>`_ (v4.1), extended to included ionizing radiative transfer, compatible with static mesh refinement and MPI parallelization.  The original Athena copyright can be found in copyright.h.  For our modification:

AUTHORS: 
Anjali Tripathi
Mark R. Krumholz

Copyright 2015 Anjali Tripathi, Mark R. Krumholz

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License v2 or later as published by the Free Software Foundation.
  
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
  
