# Satellite-Orbit-Sim

## Description
Non-linear satellite orbit simulation in Low Earth Orbit (LEO) with J2 perturbation considerations.

## Background
This was an assignment given in AA528 Space Dynamics & Control at the University of Washington by Professor Mehran Mesbahi. The problem comes from "Analytical Mechanics of Space Systems" (4th Edition) by Hanspeter Schaub and John L. Junkins. All analysis provided originates from this source.

## Problem Statement
### 12.6
Assume an Earth-relative orbit is given by the initial orbit elements _a_ = 7300 km, _e_ = 0.05, _i_ = 42 deg, &Omega; = 0 deg, &omega; = 45 deg, and _M_<sub>0</sub> = 0 deg. Assume the disturbance acceleration is solely due to the J<sub>2</sub> gravitational acceleration given in Eq. (11.57).

1. Set up a numerical simulation to solve the true nonlinear motion {_x(t), x'(t)_} using Eq. (12.1) for 10 orbits.
2. Translate the {_x(t), x'(t)_} coordinates into the corresponding classical orbit elements {_a_, _e_, _i_, &Omega;, &omega;, _M_<sub>0</sub>}.
3. Compare the numerically computed _e(t)_ motion to the analytically predicted instantaneous orbit element variations in Eq. (12.87) and the average orbit element variations in Eq. (12.88).

## References
Hanspeter Schaub, and John L. Junkins. Analytical Mechanics of Space Systems. American Institute Of Aeronautics And Astronautics, Inc, 2018.
