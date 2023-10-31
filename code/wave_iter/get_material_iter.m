function space = get_material_iter(sol, mat, space, iter_relax)
% Update the spatially dependent material parameters.
%
%    The material parameters are spatially dependent.
%    The material parameters are updated with a relaxation parameter.
%
%    Parameters:
%        sol (struct): latest solution of the full wave electromagnetic problem
%        mat (struct): description of the material parameters (permeability, permitivitty)
%        space (struct): existing description of the spatially dependent material parameters
%
%    Returns:
%        space (struct): updated description of the spatially dependent material parameters
%
%    Thomas Guillod - Dartmouth College.
%    2023 - MIT License.

% extract the solution
r_vec = sol.r_vec;
B_vec = sol.B_vec;
E_vec = sol.E_vec;

% evaluate the new material parameters
mu_new_vec = mat.mu_fct(B_vec);
eps_new_vec = mat.eps_fct(E_vec);

% evaluate the existing material parameters
mu_ref_vec = space.mu_fct(r_vec);
eps_ref_vec = space.eps_fct(r_vec);

% compute the material parameters with relaxation
mu_vec = (1-iter_relax).*mu_ref_vec+iter_relax.*mu_new_vec;
eps_vec = (1-iter_relax).*eps_ref_vec+iter_relax.*eps_new_vec;

% create spatially dependent material parameters
mu_fct = griddedInterpolant(r_vec, mu_vec, 'linear', 'nearest');
eps_fct   = griddedInterpolant(r_vec, eps_vec, 'linear', 'nearest');

% assign
space.mu_fct = mu_fct;
space.eps_fct = eps_fct;

end