# NBodyGravity
Simulation of N bodies interacting gravitationally

Run the script NBodyGravitySimulation

three figures will pop up. Figure 2 view is forced on Figure 1


## Optional pairwise gravity routine

`ThreeBodyEOMsPairwise.m` evaluates each unordered body pair once and updates
both accelerations. To try it, replace `ThreeBodyEOMs` in the `ode89` call in
`NBodyGravitySimulation.m` with `ThreeBodyEOMsPairwise`.
The original routine and default simulation remain available for comparison.

From the repository root, run `verify_in_matlab` in MATLAB or GNU Octave for
correctness checks and a routine-only benchmark. The script records its runtime.
Recorded validation passed in GNU Octave 10.3.0; native MathWorks MATLAB and the
complete graphics/`ode89` driver have not been tested by this contribution.
See [validation details](validation/pairwise/README.md).
