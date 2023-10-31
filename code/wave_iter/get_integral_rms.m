function rms = get_integral_rms(r, r_vec, var_vec)
% Compute the spatial RMS value.
%
%    Parameters:
%        r (float): radius of the cylinder geometry
%        r_vec (vector): vector with the spatial samples
%        var_vec (vector): vector with the variable values
%
%    Returns:
%        rms (float): computed spatial RMS value
%
%    Thomas Guillod - Dartmouth College.
%    2023 - MIT License.

% integrate
int = 2.*pi.*trapz(r_vec, abs(var_vec).^2.*r_vec);

% compute RMS
A = pi.*r.^2;
rms = sqrt(int./A);

end
