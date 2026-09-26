% Run with Jorge's unmodified NBodyGravity folder and this folder on the path.
% Executes the actual .m files in either MATLAB or GNU Octave.
% JSON results record the runtime; Octave results are not MATLAB results.
rng(20260926, 'twister');
cases = 0;
exact_matches = 0;
max_normalized_error = 0;
max_force_residual = 0;
scales = [1, 1; 1e-3, 1; 1e3, 6.67430e-11];
for N = [1, 2, 3, 4, 8, 16, 32]
    for scale_index = 1:size(scales,1)
        for trial = 1:12
            state = randn(6, N);
            state(1:3,:) = state(1:3,:)*scales(scale_index,1);
            masses = exp(randn(1, N));
            G = scales(scale_index,2);
            expected = ThreeBodyEOMs(0, state(:), masses, N, G);
            actual = ThreeBodyEOMsPairwise(0, state(:), masses, N, G);
            scale = max(1, max(abs(expected)));
            err = max(abs(expected-actual))/scale;
            assert(all(isfinite(actual)) && err <= 1e-12, ...
                'Pairwise derivative does not match the upstream derivative.');
            derivative = reshape(actual,6,N);
            forces = bsxfun(@times,derivative(4:6,:),masses);
            residual = norm(sum(forces,2))/max(1,sum(sqrt(sum(forces.^2,1))));
            assert(residual < 1e-12);
            max_force_residual = max(max_force_residual,residual);
            max_normalized_error = max(max_normalized_error,err);
            exact_matches = exact_matches+isequal(expected,actual);
            cases = cases+1;
        end
    end
end
s = [0 2; 0 0; 0 0; 1 -1; 2 -2; 3 -3];
expected = [s(4:6,:); 1.25 -0.75; 0 0; 0 0];
assert(max(abs(ThreeBodyEOMsPairwise(0,s(:),[3 5],2,1)-expected(:))) < 1e-14);
s = (1:6)';
assert(isequal(ThreeBodyEOMsPairwise(0,s,2,1,1),[4;5;6;0;0;0]));
fprintf('PASS: %d comparisons, %d exact matches, plus 2 analytic cases.\n',cases,exact_matches);

binary = [-.5 .5; 0 0; 0 0; 0 0; -sqrt(.5) sqrt(.5); 0 0];
triple = [-1 1 0; 0 0 2; 0 0 .2; 0 0 .1; -.3 .3 0; 0 0 0];
initials = {binary,triple};
durations = [2*pi/sqrt(2),.25];
trajectory_errors = zeros(1,2);
orbit_closure_error = NaN;
options = odeset('RelTol',1e-11,'AbsTol',1e-13);
for k = 1:2
    initial = initials{k};
    N = size(initial,2);
    masses = ones(1,N);
    times = linspace(0,durations(k),101);
    [~,reference] = ode45(@(t,y)ThreeBodyEOMs(t,y,masses,N,1),times,initial(:),options);
    [~,candidate] = ode45(@(t,y)ThreeBodyEOMsPairwise(t,y,masses,N,1),times,initial(:),options);
    trajectory_errors(k) = max(abs(reference(:)-candidate(:)));
    assert(trajectory_errors(k) < 1e-10);
    if k == 1
        orbit_closure_error = max(abs(candidate(end,:)'-initial(:)));
        assert(orbit_closure_error < 1e-8);
    end
end
fprintf('PASS: 2 ode45 trajectories; differences %.3g, %.3g.\n',trajectory_errors);

% Kernel-only benchmark: warm-up, then seven alternating-order batches.
% Repetition count targets >=0.10 seconds for the original routine per batch.
benchmark_sizes = [3,16,64,128];
benchmarks = struct([]);
for k = 1:length(benchmark_sizes)
    N = benchmark_sizes(k);
    state = randn(6,N);
    masses = exp(randn(1,N));
    old_fn = @()ThreeBodyEOMs(0,state(:),masses,N,1);
    new_fn = @()ThreeBodyEOMsPairwise(0,state(:),masses,N,1);
    old_fn(); new_fn();
    repetitions = 1;
    while true
        tic;
        for j = 1:repetitions, old_fn(); end
        elapsed = toc;
        if elapsed >= .10 || repetitions >= 4096, break; end
        repetitions = repetitions*2;
    end
    old_times = zeros(1,7); new_times = zeros(1,7);
    for batch = 1:7
        functions = {old_fn,new_fn};
        if mod(batch,2) == 0, order = [2,1]; else, order = [1,2]; end
        for which = order
            fn = functions{which};
            tic;
            for j = 1:repetitions, fn(); end
            seconds_per_call = toc/repetitions;
            if which == 1, old_times(batch) = seconds_per_call;
            else, new_times(batch) = seconds_per_call; end
        end
    end
    benchmarks(k).bodies = N;
    benchmarks(k).repetitions_per_batch = repetitions;
    benchmarks(k).original_seconds = old_times;
    benchmarks(k).candidate_seconds = new_times;
    benchmarks(k).original_median_seconds = median(old_times);
    benchmarks(k).candidate_median_seconds = median(new_times);
    benchmarks(k).median_time_ratio = median(old_times)/median(new_times);
    fprintf('N=%d: original %.6g s; pair-once %.6g s; ratio %.3gx\n', ...
        N,median(old_times),median(new_times),benchmarks(k).median_time_ratio);
end
is_octave = exist('OCTAVE_VERSION','builtin') ~= 0;
if is_octave, runtime_name = 'GNU Octave'; else, runtime_name = 'MathWorks MATLAB'; end
report = struct('status','PASS','runtime',runtime_name,'runtime_version',version, ...
    'native_matlab_executed',~is_octave,'seed',20260926, ...
    'upstream_comparison_cases',cases,'exact_numeric_matches',exact_matches, ...
    'max_normalized_derivative_difference',sprintf('%.17g',max_normalized_error), ...
    'max_normalized_total_force_residual',sprintf('%.17g',max_force_residual), ...
    'analytic_cases_passed',2,'trajectory_solver','ode45', ...
    'trajectory_max_absolute_differences',trajectory_errors, ...
    'circular_binary_closure_error',orbit_closure_error,'benchmarks',benchmarks);
if is_octave, filename = 'octave_validation_results.json';
else, filename = 'matlab_validation_results.json'; end
fid = fopen(fullfile(fileparts(mfilename('fullpath')),filename),'w');
assert(fid >= 0);
fprintf(fid,'%s\n',jsonencode(report));
fclose(fid);
fprintf('Saved %s (%s %s).\n',filename,runtime_name,version);
