function x_noisy = add_noise_state(x, N1, N2, noise_level)
%ADD_NOISE_STATE  Apply multiplicative Gaussian noise to the epidemiological state.
%
% This function perturbs the 6-dimensional state vector
%   x = [S1; V1; I1; S2; V2; I2]
% by a multiplicative random factor on each component:
%   x_noisy = x .* scale,
% where scale_i = 1 + noise_level * Z_i,   Z_i ~ N(0,1).
%
%noise_level is the standard deviation of the multiplicative perturbation.
%
% Renormalization:
%  - After perturbation, the function enforces that the totals in each age class do not
%    exceed the corresponding population sizes N1 and N2:
%       S1+V1+I1 <= N1,    S2+V2+I2 <= N2.
%  - If a total exceeds N1 (resp. N2), the three components are rescaled proportionally.
%  - If a total is below N1 (resp. N2), the "missing mass" is implicitly assigned to the
%    recovered compartment R1 (resp. R2), which is not stored explicitly.
%
% Inputs:
%  - x          : 6x1 state vector [S1;V1;I1;S2;V2;I2]
%  - N1, N2     : population sizes of infants and non-infants
%  - noise_level: standard deviation of multiplicative noise (e.g., 0.05)
%
% Output:
%  - x_noisy    : perturbed and (if needed) renormalized state vector

x = x(:);

% Draw a 6x1 vector of independent standard Gaussian variables.
% The multiplicative scaling factors have mean 1 and standard deviation noise_level:
%   scale_i = 1 + noise_level * Z_i,   Z_i ~ N(0,1).
% Roughly, by the "3-sigma rule", scale_i will lie in [1-3*noise_level, 1+3*noise_level]
% with high probability (but not guaranteed).
scale = 1 + noise_level * randn(6,1);

% Clamp scale to the interval [1-3*noise_level, 1+3*noise_level] to avoid extreme values.
scale = max(1-3*noise_level, min(1+3*noise_level, scale));

% Apply multiplicative perturbation component-wise
x_noisy = x .* scale;

%% Renormalization by age class
% Age class 1 (infants): indices 1:3 correspond to [S1; V1; I1]
tot1 = sum(x_noisy(1:3));
if tot1 > N1
    % Rescale proportionally so that S1+V1+I1 equals N1
    x_noisy(1:3) = x_noisy(1:3) * (N1 / tot1);
end
% If tot1 < N1, the remainder is implicitly in R1 (not stored explicitly)

% Age class 2 (non-infants): indices 4:6 correspond to [S2; V2; I2]
tot2 = sum(x_noisy(4:6));
if tot2 > N2
    % Rescale proportionally so that S2+V2+I2 equals N2
    x_noisy(4:6) = x_noisy(4:6) * (N2 / tot2);
end
% If tot2 < N2, the remainder is implicitly in R2

end