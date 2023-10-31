function [sol, is_valid] = get_solve_approx(problem, mat)
% Solve the magnetic-dielectric effects in a core with a quasi-static approximation.
%
%    A infinite cylindrical core is considered.
%    A quasi-static approximation is used (no wave propagation).
%    The (average) flux density through the cross-section is imposed (excitation).
%    The material parameters can be amplitude-dependent (locally linearized).
%    The material parameters are constant over the cross-section (spatially independent).
%
%    Parameters:
%        problem (struct): description of the problem (frequency, flux density, radius)
%        mat (struct): description of the material parameters (permeability, permitivitty)
%
%    Returns:
%        sol (struct): solution of the quasi-static electromagnetic problem
%        is_valid (boolean): validity of the solution (detect extrapolation)
%
%    Thomas Guillod - Dartmouth College.
%    2023 - MIT License.

% extract problem data
f = problem.f;
r = problem.r;
B_src = problem.B_src;

% compute the quasi-static approxiation
B = B_src;
E = (2.*pi.*f.*pi.*r.^2.*B_src)./(2.*pi.*r);
E = E./sqrt(2);

% evaluate the material parameters
[mu, is_valid_mu] = mat.mu_fct(B);
[eps, is_valid_eps] = mat.eps_fct(E);

% check if the material parameters are extrapolated
is_valid = is_valid_mu&is_valid_eps;

% evaluate the fields
[H, D, J, dBdt] = get_post_field(f, B, E, mu, eps);

% compute the complex power
[s_mag, s_ele, s_tot] = get_post_power(E, J, H, dBdt);

% compute the wavelength and penetration depth
[lambda, delta] = get_post_wavenumber(f, mu, eps);

% assign
sol.mu = mu;
sol.eps = eps;
sol.E = abs(E);
sol.B = abs(B);
sol.H = abs(H);
sol.D = abs(D);
sol.J = abs(J);
sol.s_mag = s_mag;
sol.s_ele = s_ele;
sol.s_tot = s_tot;
sol.lambda = lambda;
sol.delta = delta;

end
