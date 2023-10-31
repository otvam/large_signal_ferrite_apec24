function avg = get_integral_avg(r, r_vec, var_vec)
% Compute the spatial average value.
%
%    Parameters:
%        r (float): radius of the cylinder geometry
%        r_vec (vector): vector with the spatial samples
%        var_vec (vector): vector with the variable values
%
%    Returns:
%        avg (float): computed spatial average value
%
%    Thomas Guillod - Dartmouth College.
%    2023 - MIT License.

% integrate
int = 2.*pi.*trapz(r_vec, var_vec.*r_vec);

% compute average
A = pi.*r.^2;
avg = int./A;

end
