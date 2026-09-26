function ds = ThreeBodyEOMsPairwise(t, s, m, N, G)
% Pair-once Newtonian gravity RHS; same state layout as ThreeBodyEOMs.
% Contribution from Ask Badger (jdbjdbc-badger), with ChatGPT assistance.
% Original project: jorge-physics/NBodyGravity.
% Intended domain: finite real double inputs, N >= 1, distinct positions.
% No approximation, force cutoff, softening, or integrator change is added.
% Direct GNU Octave validation passed; MathWorks MATLAB was not run.
% See verify_in_matlab.m and the recorded validation results.

state = reshape(s, 6, N);
acceleration = zeros(3, N);
for i = 1:N-1
    for j = i+1:N
        displacement = state(1:3,j) - state(1:3,i);
        distance_cubed = norm(displacement)^3;
        acceleration(:,i) = acceleration(:,i) + ...
            (G*m(j)/distance_cubed).*displacement;
        acceleration(:,j) = acceleration(:,j) - ...
            (G*m(i)/distance_cubed).*displacement;
    end
end
derivative = [state(4:6,:); acceleration];
ds = derivative(:);
end
