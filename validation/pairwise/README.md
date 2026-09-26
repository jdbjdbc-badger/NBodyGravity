# Pairwise gravity validation

Contribution from Ask Badger / jdbjdbc-badger, prepared with ChatGPT assistance.
Pair symmetry is a standard Newtonian identity applied to this implementation.
No private Ask Badger engine code is included.

The candidate calculates each unordered pair's distance once, updates both
bodies, removes the N-by-N coefficient table, and avoids repeated setdiff calls.
It preserves all interactions and the existing six-values-per-body state layout.
There is no cutoff, softening, approximation, or integrator change.
Work remains O(N^2); auxiliary working storage becomes O(N).

## Recorded checks

The retained September 26, 2026 run executed the actual original and candidate
MATLAB-language files in GNU Octave 10.3.0 on Linux x86-64:

- 252 randomized derivative comparisons, all numerically identical.
- Two analytical cases: unequal-mass binary and free single body.
- Two ode45 trajectory comparisons, identical at all 101 sampled times.
- Circular-binary one-orbit closure error: 1.9446e-11 in test units.
- Maximum normalized net-force residual: 1.2334843288814213e-16.

Seven alternating-order timing batches followed warm-up. Median original/candidate
routine-time ratios were 20.20x (3 bodies), 5.83x (16), 2.78x (64), and 2.44x (128).
Raw timings and the run log are included here. These are Octave routine timings,
not native MATLAB or whole-application speedups. All organization changes contribute
to the result; the ratios do not isolate pair reuse alone.

Test domain: finite real doubles, N >= 1, distinct positions, no overflow.
Singular collisions and long chaotic trajectories were not validated.
MathWorks MATLAB was not available. The trajectory tests used ode45, not the
driver's ode89 solver, and did not exercise graphics or CSV output.

## Reproduce

From the repository root in MATLAB or Octave:

```matlab
verify_in_matlab
```

The script writes a runtime-labelled JSON file beside itself. Results in this
directory are the retained Octave run, not results from the reader's machine.
The original source identity is recorded in upstream_source.json. The published
candidate differs from the tested candidate only in its attribution comment.
