function [r_vec, v_vec, d_vec] = get_helmholtz_wave(r, f, bnd, param, space)
% Solve the Helmholtz wave equation.
%
%    Parameters:
%        r (float): radius of the cylinder geometry
%        f (float): frequency for solving the equation
%        bnd (float): value of the boundary condition
%        param (struct): description of the solver numerical parameters and tolerances
%        space (struct): description of the spatially dependent material parameters
%
%    Returns:
%        r_vec (vector): vector describing the sampling along the cross-section
%        v_vec (vector): vector with the solution (cylindrical coordinate)
%        d_vec (vector): vector with the derivative (cylindrical coordinate)
%
%    Thomas Guillod - Dartmouth College.
%    2023 - MIT License.

% extract solver parameters
tol_pde = param.tol_pde;
tol_int = param.tol_int;
n_init = param.n_init;
n_sample = param.n_sample;

% function describing the PDE
fct_pde = @(r_vec, y_vec) get_pde(r_vec, y_vec, f, space);

% function describing the BC
fct_bnd = @(y_a, y_b) get_bnd(y_a, y_b, bnd);

% function returning the initial guess
fct_init = @(r_vec) get_init(r_vec, r, bnd);

% vector describing the sampling along the cross-section
r_vec = linspace(0, r, n_init);

% get the initial guess
sol_init = bvpinit(r_vec, fct_init);

% options for the PDE solver
options = bvpset('RelTol', tol_pde,'Stats','on', 'Vectorized', 'on', 'Stats', 'off');

% solve the PDE
sol = bvp4c(fct_pde, fct_bnd, sol_init, options);

% extract the solution
r_sol = sol.x;
v_sol = sol.y(1,:);
dv_sol = sol.y(2,:);

% interpolant for the solution
fct_v = griddedInterpolant(r_sol, v_sol, 'linear', 'linear');
fct_dv = griddedInterpolant(r_sol, dv_sol, 'linear', 'linear');

% resample the solution wth the desired samples
r_vec = linspace(0, r, n_sample);
v_vec = fct_v(r_vec);
dv_vec = fct_dv(r_vec);

% get solution derivative in polar coordinate
d_vec_singular = 2.*dv_vec;
d_vec_normal = dv_vec+v_vec./r_vec;

% fix the singularity of the derivative at the origin
idx = r_vec<eps;
d_vec(idx==true) = d_vec_singular(idx==true);
d_vec(idx==false) = d_vec_normal(idx==false);

% check that the integral of the derivative matches the boundary condition
int_ana = bnd.*2.*pi.*r;
int_fem = 2.*pi.*abs(trapz(r_vec, r_vec.*d_vec));
err_int = abs(int_fem-int_ana)./int_ana;
assert(err_int<tol_int, 'solver convergence issue')

end

function dy_vec = get_pde(r_vec, y_vec, f, space)
% Function describing the PDE.
%
%    Parameters:
%        r_vec (vector): vector with the spatial samples
%        y_vec (vector): vector with the solution
%        f (float): frequency for solving the equation
%        space (struct): description of the spatially dependent material parameters
%
%    Returns:
%        dy_vec (vector): vector with the solution derivative

% extract the solution
v_vec = y_vec(1, :);
dv_vec = y_vec(2, :);

% constant
mu_0 = 4.*pi.*1e-7;
eps_0 = 8.854.*1e-12;

% evaluate the spatially dependent material parameters
mu_vec = mu_0.*space.mu_fct(r_vec);
eps_vec = eps_0.*space.eps_fct(r_vec);

% compute the wavenumbers
w = 2.*pi.*f;
k_vec = w.*sqrt(mu_vec.*eps_vec);

% compute the square of the samples and wavenumbers
r2_vec = r_vec.^2;
k2_vec = k_vec.^2;

% derivative for the non-singular points
dy_normal_vec = [dv_vec ; -(1./r_vec).*dv_vec+(1./r2_vec).*v_vec-k2_vec.*v_vec];

% derivative for the singular points (at the origin)
dy_singular_vec = [dv_vec ; -k2_vec.*v_vec];

% assemble the derivative
idx = r_vec<eps;
dy_vec(:, idx==true) = dy_singular_vec(:, idx==true);
dy_vec(:, idx==false) = dy_normal_vec(:, idx==false);

end

function res = get_bnd(y_a, y_b, bnd)
% Function describing the BC.
%
%    Parameters:
%        y_a (vector): solution at the origin
%        y_b (vector): solution at the boundary
%        bnd (float): value of the boundary condition
%
%    Returns:
%        res (vector): residuum of the boundary condition

% extract
v_a = y_a(1, :);
v_b = y_b(1, :);

% compute the residuum
res = [v_a ; v_b-bnd];

end

function y_vec = get_init(r_vec, r, bnd)
% Function returning the initial guess.
%
%    Parameters:
%        r_vec (vector): vector with the spatial samples
%        r (float): radius of the cylinder geometry
%        bnd (float): value of the boundary condition
%
%    Returns:
%        y_vec (vector): vector with the initial guess

% get initial guess (low-frequency approximation)
shape = ones(1, length(r_vec));
v_vec = (r_vec.*bnd./r).*shape;
dv_vec = (bnd./r).*shape;

% assemble the solution
y_vec = [v_vec ; dv_vec];

end