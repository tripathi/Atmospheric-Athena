Atmospheric Athena:
==================
3D Atmospheric escape model with ionizing radiative transfer

What is Atmospheric Athena?
  Atmospheric Athena is a code intended to simulate hydrodynamic escape from close-in giant planets in 3D.  It uses the Athena hydrodynamics code v4.1 (Stone et al. 2008) with a new ionizing radiative transfer implementation based on Krumholz et al, 2007, to self-consistently model photoionization driven winds from the planet.  The code is fully compatible with static mesh refinement and MPI parallelization.  It can handle arbitrary planet potentials and stellar initial conditions.  The physical system that the code models and the equations that it solves are fully described in Tripathi, Kratter, Murray-Clay, & Krumholz, submitted 2015.