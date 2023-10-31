function [sol, is_valid] = get_solve_single(problem, param, mat, space)
% Solve the electromagnetic problem with given spatially dependent material parameters.
%
%    Parameters:
%        problem (struct): description of the problem (frequency, flux density, radius)
%        param (struct): description of the solver numerical parameters and tolerances
%        mat (struct): description of the material parameters (permeability, permitivitty)
%        space (struct): existing description of the spatially dependent material parameters
%
%    Returns:
%        sol (struct): solution of the full wave electromagnetic problem
%        is_valid (boolean): validity of the solution (detect extrapolation)
%
%    Thomas Guillod - Dartmouth College.
%    2023 - MIT License.

% extract problem data
f = problem.f;
r = problem.r;
B_src = problem.B_src;

% extract solver parameters
frac_valid_eps = param.frac_valid_eps;
frac_valid_mu = param.frac_valid_mu;

% find the boundary condition
bnd = (2.*pi.*f.*pi.*r.^2.*B_src)./(2.*pi.*r);

% solve the Helmholtz wave equation
[r_vec, v_vec, d_vec] = get_helmholtz_wave(r, f, bnd, param, space);

% assign the solution of the wave equation
E_vec = v_vec;
B_vec = 1i.*d_vec./(2.*pi.*f);

% evaluate the material parameters
[mu_vec, is_valid_mu_vec] = mat.mu_fct(B_vec);
[eps_vec, is_valid_eps_vec] = mat.eps_fct(E_vec);

% check if the material parameters are extrapolated
frac_mu = nnz(is_valid_mu_vec)./numel(is_valid_mu_vec);
frac_eps= nnz(is_valid_eps_vec)./numel(is_valid_eps_vec);
is_valid = (frac_mu>=frac_valid_mu)&(frac_eps>=frac_valid_eps);

% evaluate the fields
[H_vec, D_vec, J_vec, dBdt_vec] = get_post_field(f, B_vec, E_vec, mu_vec, eps_vec);

% compute the complex power
[s_mag_vec, s_ele_vec, s_tot_vec] = get_post_power(E_vec, J_vec, H_vec, dBdt_vec);

% compute the wavelength and penetration depth
[lambda_vec, delta_vec] = get_post_wavenumber(f, mu_vec, eps_vec);

% compute the spatial RMS value of the fields
B = get_integral_rms(r, r_vec, B_vec);
E = get_integral_rms(r, r_vec, E_vec);
H = get_integral_rms(r, r_vec, H_vec);
D = get_integral_rms(r, r_vec, D_vec);
J = get_integral_rms(r, r_vec, J_vec);

% compute the average complex power
s_mag = get_integral_avg(r, r_vec, s_mag_vec);
s_ele = get_integral_avg(r, r_vec, s_ele_vec);
s_tot = get_integral_avg(r, r_vec, s_tot_vec);

% compute the average material parameters
mu = get_integral_avg(r, r_vec, mu_vec);
eps = get_integral_avg(r, r_vec, eps_vec);

% compute the average wavelength and penetration depth
lambda = get_integral_avg(r, r_vec, lambda_vec);
delta = get_integral_avg(r, r_vec, delta_vec);

% assign scalar
sol.mu = mu;
sol.eps = eps;
sol.E = E;
sol.B = B;
sol.H = H;
sol.D = D;
sol.J = J;
sol.s_mag = s_mag;
sol.s_ele = s_ele;
sol.s_tot = s_tot;
sol.lambda = lambda;
sol.delta = delta;

% assign vector
sol.r_vec = r_vec;
sol.mu_vec = mu_vec;
sol.eps_vec = eps_vec;
sol.E_vec = E_vec;
sol.B_vec = B_vec;
sol.H_vec = H_vec;
sol.D_vec = D_vec;
sol.J_vec = J_vec;
sol.s_mag_vec = s_mag_vec;
sol.s_ele_vec = s_ele_vec;
sol.lambda_vec = lambda_vec;
sol.delta_vec = delta_vec;

end

