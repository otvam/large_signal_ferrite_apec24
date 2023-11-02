function [problem, param, mat, s_ref] = get_param(f, B_src, A)
% Get the solver parameters for a given electromagnetic problem.
%
%    Parameters:
%        f (float): frequency
%        B_src (float): average flux density in the core
%        A (float): cross-section of the magnetic core
%
%    Returns:
%        problem (struct): description of the problem (frequency, flux density, radius)
%        param (struct): description of the solver numerical parameters and tolerances
%        mat (struct): description of the material parameters (permeability, permitivitty)
%        s_ref (float): complex power without dielectric effects
%
%    Thomas Guillod - Dartmouth College.
%    2023 - MIT License.


% constant
mu_0 = 4.*pi.*1e-7;

% description of the problem (frequency, flux density, radius)
problem.f = f;              % frequency
problem.B_src = B_src;      % average flux density in the core
problem.r = sqrt(A./pi);    % radius of the magnetic core

% description of the solver numerical parameters and tolerances
param.tol_pde = 1e-5;          % relative tolerance for the PDE solver
param.tol_int = 1e-5;          % relative tolerance for checking the BC
param.n_init = 500;            % number of initial samples for the PDE solver
param.n_sample = 5000;         % number of samples for the solution
param.frac_valid_mu = 0.8;     % required fraction of non-extrapolated points for the permeability
param.frac_valid_eps = 0.8;    % required fraction of non-extrapolated points for the permitivitty
param.iter_min = 3;            % minimum number of solver iteration
param.iter_tol = 0.2e-2;       % iteration relative tolerance (on the complex power)
param.iter_relax = 0.1;        % iteration relaxation parameter (for the material parameters)

% description of the material parameters (permeability, permitivitty)
mu_fct = get_mu_fct('large_signal', 'N87_mu_large', f);       % large-signal complex permeability
eps_fct = get_eps_fct('large_signal', 'N87_eps_large', f);    % large-signal complex permitivitty
mat = struct('mu_fct', mu_fct, 'eps_fct', eps_fct);

% compute the complex power without dielectric effects
mu = mu_fct(B_src);
s_ref = 0.5.*2.*pi.*f.*(1i./conj(mu_0.*mu)).*(B_src.^2);

end