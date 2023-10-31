function [lambda, delta] = get_post_wavenumber(f, mu, eps)
% Compute the wavelength and penetration depth.
%
%    Parameters:
%        f (float): frequency
%        mu (float): relative permeability
%        eps (float): relative permitivitty
%
%    Returns:
%        lambda (float): wavelength
%        delta (float): penetration depth
%
%    Thomas Guillod - Dartmouth College.
%    2023 - MIT License.

% constant
mu_0 = 4.*pi.*1e-7;
eps_0 = 8.854e-12;

% angular frequency
w = 2.*pi.*f;

% compute wavenumber
k_square = w.*mu_0.*mu.*w.*eps_0.*eps;
k = sqrt(k_square);

% split real and imaginary
beta = +real(k);
alpha = -imag(k);

% assign
lambda = 2.*pi./beta;
delta = 1./alpha;

end