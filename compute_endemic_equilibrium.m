function [x_star, info] = compute_endemic_equilibrium(par, x_seed)
%COMPUTE_ENDEMIC_EQUILIBRIUM Find and diagnose the positive fixed point.
%   The weekly model is written as x(t+1) = F(x(t)). This function solves
%   F(x*) - x* = 0 after a long burn-in that provides an initial guess in
%   the basin of attraction of the endemic, rather than disease-free,
%   equilibrium. A log transformation keeps all six explicit states
%   positive during the nonlinear solve.

    if numel(x_seed) ~= 6 || any(~isfinite(x_seed)) || any(x_seed <= 0)
        error('compute_endemic_equilibrium: x_seed must contain six positive finite states.');
    end

    burnin_weeks = 5000;
    [x_burnin, ~, ~, ~, ~, ~, ~, ~, ~, ~] = ...
        discrete_system(0, burnin_weeks, x_seed, par);
    x_guess = x_burnin(:,end);

    % Solve in log coordinates: x = x_scale .* exp(z), with z0 = 0.
    x_scale = max(x_guess, 1);
    residual = @(z) local_scaled_log_residual(z, x_scale, par);
    options = optimoptions('fsolve', ...
        'Display', 'off', ...
        'FunctionTolerance', 1e-12, ...
        'StepTolerance', 1e-12, ...
        'OptimalityTolerance', 1e-12, ...
        'MaxIterations', 1000, ...
        'MaxFunctionEvaluations', 10000);

    [z_star, ~, exitflag, output] = fsolve(residual, zeros(6,1), options);
    x_star = x_scale .* exp(z_star);

    [x_next, ~, ~, inc1, inc2, ~, ~, ~, R1, R2] = ...
        discrete_system(0, 1, x_star, par);
    fixed_point_residual = x_next(:,end) - x_star;

    R1_star = R1(1);
    R2_star = R2(1);
    if exitflag <= 0
        error('compute_endemic_equilibrium: fsolve failed (exitflag %d).', exitflag);
    end
    if any(x_star <= 0) || R1_star < -1e-6 || R2_star < -1e-6
        error('compute_endemic_equilibrium: the computed fixed point is not epidemiologically admissible.');
    end
    if x_star(3) < 1e-8 || x_star(6) < 1e-8
        warning('compute_endemic_equilibrium: the solver converged to a disease-free fixed point.');
    end

    map_jacobian = local_map_jacobian(x_star, par);
    jacobian_eigenvalues = eig(map_jacobian);
    spectral_radius = max(abs(jacobian_eigenvalues));

    info = struct( ...
        'R1', R1_star, ...
        'R2', R2_star, ...
        'weekly_incidence', [inc1(1); inc2(1)], ...
        'annual_incidence', 52 * [inc1(1); inc2(1)], ...
        'fixed_point_residual_inf', norm(fixed_point_residual, Inf), ...
        'fixed_point_residual_relative', max(abs(fixed_point_residual) ./ max(1, abs(x_star))), ...
        'map_jacobian', map_jacobian, ...
        'jacobian_eigenvalues', jacobian_eigenvalues, ...
        'spectral_radius', spectral_radius, ...
        'locally_stable', spectral_radius < 1, ...
        'burnin_weeks', burnin_weeks, ...
        'exitflag', exitflag, ...
        'solver_output', output);
end

function res = local_scaled_log_residual(z, x_scale, par)
    x = x_scale .* exp(z);
    [x_next, ~, ~, ~, ~, ~, ~, ~, ~, ~] = discrete_system(0, 1, x, par);
    res = (x_next(:,end) - x) ./ x_scale;
end

function J = local_map_jacobian(x, par)
    n = numel(x);
    J = zeros(n);

    for j = 1:n
        h = 1e-6 * max(1, abs(x(j)));
        x_plus = x;
        x_minus = x;
        x_plus(j) = x_plus(j) + h;
        x_minus(j) = x_minus(j) - h;

        F_plus = local_one_step_map(x_plus, par);
        if x_minus(j) > 0
            F_minus = local_one_step_map(x_minus, par);
            J(:,j) = (F_plus - F_minus) / (2*h);
        else
            F_base = local_one_step_map(x, par);
            J(:,j) = (F_plus - F_base) / h;
        end
    end
end

function F = local_one_step_map(x, par)
    [x_next, ~, ~, ~, ~, ~, ~, ~, ~, ~] = discrete_system(0, 1, x, par);
    F = x_next(:,end);
end
