function space = get_material_init(problem, mat)
% Init the spatially dependent material parameters for the first iteration.
%
%    For the first iteration, the material parameters are constant.
%    For the first iteration, the quasi-static approximation is used.
%
%    Parameters:
%        problem (struct): description of the problem (frequency, flux density, radius)
%        mat (struct): description of the material parameters (permeability, permitivitty)
%
%    Returns:
%        space (struct): created description of the spatially dependent material parameters
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
mu = mat.mu_fct(abs(B));
eps = mat.eps_fct(abs(E));

% create spatially independent material parameters
mu_fct = @(r) mu.*ones(size(r));
eps_fct = @(r) eps.*ones(size(r));

% assign
space.mu_fct = mu_fct;
space.eps_fct = eps_fct;

end
