function mu_out_fct = get_mu_fct(type, dataset, f)
% Return a function describing the material complex permeability for a given frequency.
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
%        mu_out_fct (function): function describing the amplitude dependent material parameters
%
%    Thomas Guillod - Dartmouth College.
%    2023 - MIT License.

% load dataset
data = load(['dataset/' dataset '/data.mat']);

% get the function
%     - describing the complex permeability
%     - describing the extrapolation
switch type
    case 'small_signal'
        % load data
        f_vec = data.f_vec;
        mu_vec = data.mu_vec;

        % interpolate frequency
        mu_tmp = interp1(log10(f_vec), mu_vec, log10(f));

        % check validity
        assert(isfinite(mu_tmp), 'invalid frequency')

        % get the interpolant (independent of the amplitude)
        mu_fct = @(B) mu_tmp.*ones(size(B));
        is_valid_fct = @(B) true(size(B));
    case 'large_signal'
        % load data
        f_vec = data.f_vec;
        B_mat = data.B_mat;
        mu_mat = data.mu_mat;

        % interpolate frequency
        B_vec = 10.^interp1(log10(f_vec), log10(B_mat), log10(f));
        mu_vec = interp1(log10(f_vec), mu_mat, log10(f));

        % get the interpolant (dependent of the amplitude)
        [mu_fct, is_valid_fct] = get_interp(B_vec, mu_vec);
    otherwise
        error('invalid tag')
end

% assemble the functions into a single function handle
mu_out_fct = @(B) get_out(B, mu_fct, is_valid_fct);

end

function [mu_vec, is_valid_vec] = get_out(B_vec, mu_fct, is_valid_fct)
% Function returning the complex material parameters and the extrapolation flag.
%
%    Parameters:
%        B_vec (vector): vector with the flux density amplitude
%        mu_fct (function): function describing the complex permeability
%        is_valid_fct (function): function describing extrapolation
%
%    Returns:
%        mu_vec (vector): evaluated complex permeability
%        is_valid_vec (vector): flag indicating if the material parameters are extrapolated 

mu_vec = mu_fct(B_vec);
is_valid_vec = is_valid_fct(B_vec);

end

function [mu_fct, is_valid_fct] = get_interp(B_vec, mu_vec)
% Create interpolation for the material parameters and the extrapolation flag.
%
%    Parameters:
%        B_vec (vector): vector with the measured flux density points
%        mu_vec (vector): vector with the measured complex permeability
%
%    Returns:
%        mu_fct (function): function describing the complex permeability
%        is_valid_fct (function): function describing extrapolation

% filter invalid points
idx = isfinite(B_vec)&isfinite(mu_vec);
B_vec = B_vec(idx);
mu_vec = mu_vec(idx);

% check that interpolation is possible
assert(nnz(idx)>=2, 'invalid frequency')

% detect extrapolation
is_valid_fct = @(B) (abs(B)>=min(B_vec))&(abs(B)<=max(B_vec));

% create the interpolant (clamp for extrapolation)
mu_tmp = griddedInterpolant(log10(B_vec), mu_vec, 'linear', 'nearest');
mu_fct = @(B) mu_tmp(log10(abs(B)));

end
