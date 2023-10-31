function eps_out_fct = get_eps_fct(type, dataset, f)
% Return a function describing the material complex permitivitty for a given frequency.
%
%    With small-signal measurements, the material parameters are independent of the amplitude.
%    With large-signal measurements, the material parameters are dependent of the amplitude.
%    The function return a boolean flag indicating if the material parameters are extrapolated.
%    
%    Parameters:
%        type (str): type of the material dataset to be used (small-signal or large-signal)
%        dataset (str): name of the material dataset to be used (material name)
%        f (float): frequency used for evaluating the material parameters
%
%    Returns:
%        eps_out_fct (function): function describing the amplitude dependent material parameters
%
%    Thomas Guillod - Dartmouth College.
%    2023 - MIT License.

% load dataset
data = load(['dataset/' dataset '/data.mat']);

% get the function
%     - describing the real permitivitty
%     - describing the real conductivity
%     - describing the extrapolation
switch type
    case 'small_signal'
        % load data
        f_vec = data.f_vec;
        eps_vec = data.eps_vec;
        sigma_vec = data.sigma_vec;

        % interpolate frequency
        eps = interp1(log10(f_vec), eps_vec, log10(f));
        sigma = interp1(log10(f_vec), sigma_vec, log10(f));

        % check validity
        assert(isfinite(eps), 'invalid frequency')
        assert(isfinite(sigma), 'invalid frequency')

        % get the interpolant (independent of the amplitude)
        eps_fct = @(E) eps.*ones(size(E));
        sigma_fct = @(E) sigma.*ones(size(E));
        is_valid_fct = @(E) true(size(E));
    case 'large_signal'
        % load data
        f_vec = data.f_vec;
        E_mat = data.E_mat;
        eps_mat = data.eps_mat;
        sigma_mat = data.sigma_mat;

        % interpolate frequency
        E_vec = interp1(log10(f_vec), E_mat, log10(f));
        eps_vec = interp1(log10(f_vec), eps_mat, log10(f));
        sigma_vec = interp1(log10(f_vec), sigma_mat, log10(f));

        % get the interpolant (dependent of the amplitude)
        [eps_fct, sigma_fct, is_valid_fct] = get_interp(E_vec, eps_vec, sigma_vec);
    otherwise
        error('invalid tag')
end

% assemble the functions into a single function handle
eps_out_fct = @(E) get_out(E, eps_fct, sigma_fct, is_valid_fct, f);

end

function [eps_vec, is_valid_vec] = get_out(E_vec, eps_fct, sigma_fct, is_valid_fct, f)
% Function returning the complex material parameters and the extrapolation flag.
%
%    Parameters:
%        E_vec (vector): vector with the electric field amplitude
%        eps_fct (function): function describing the real permitivitty
%        sigma_fct (function): function describing the real conductivity
%        is_valid_fct (function): function describing extrapolation
%
%    Returns:
%        eps_vec (vector): evaluated complex permitivitty
%        is_valid_vec (vector): flag indicating if the material parameters are extrapolated 

% constant
eps_0 = 8.854.*1e-12;

% evaluate the interpolant
eps_vec = eps_fct(E_vec);
sigma_vec = sigma_fct(E_vec);
is_valid_vec = is_valid_fct(E_vec);

% transform the real permitivitty and conductivity into the complex permitivitty
eps_vec = eps_vec-1i.*sigma_vec./(2.*pi.*f.*eps_0);

end

function [eps_fct, sigma_fct, is_valid_fct] = get_interp(E_vec, eps_vec, sigma_vec)
% Create interpolation for the material parameters and the extrapolation flag.
%
%    Parameters:
%        E_vec (vector): vector with the measured electric field points
%        eps_vec (vector): vector with the measured real permitivitty
%        sigma_vec (vector): vector with the measured real conductivity
%
%    Returns:
%        eps_fct (function): function describing the real permitivitty
%        sigma_fct (function): function describing the real conductivity
%        is_valid_fct (function): function describing extrapolation

% filter invalid points
idx = isfinite(E_vec)&isfinite(eps_vec)&isfinite(sigma_vec);
E_vec = E_vec(idx);
eps_vec = eps_vec(idx);
sigma_vec = sigma_vec(idx);

% check that interpolation is possible
assert(nnz(idx)>=2, 'invalid frequency')

% detect extrapolation
is_valid_fct = @(E) (abs(E)>=min(E_vec))&(abs(E)<=max(E_vec));

% create the interpolant (clamp for extrapolation)
eps_fct_tmp = griddedInterpolant(E_vec, eps_vec, 'linear', 'nearest');
sigma_fct_tmp = griddedInterpolant(E_vec, sigma_vec, 'linear', 'nearest');
eps_fct = @(E) eps_fct_tmp(abs(E));
sigma_fct = @(E) sigma_fct_tmp(abs(E));

end
