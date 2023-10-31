function [H, D, J, dBdt] = get_post_field(f, B, E, mu, eps)
% Compute the electromagnetic field values with the material parameters.
%
%    Parameters:
%        f (float): frequency
%        B (float): magnetic flux density
%        E (float): electric field
%        mu (float): relative permeability
%        eps (float): relative permitivitty
%
%    Returns:
%        H (float): magnetic field
%        D (float): displacement field
%        J (float): current density
%        dBdt (float): magnetic flux density derivative
%
%    Thomas Guillod - Dartmouth College.
%    2023 - MIT License.

% constant
mu_0 = 4.*pi.*1e-7;
eps_0 = 8.854.*1e-12;

% apply the material relations
H = B./(mu_0.*mu);
D = eps_0.*eps.*E;

% get the field derivatives
J = 1i.*2.*pi.*f.*D;
dBdt = 1i.*2.*pi.*f.*B;

end