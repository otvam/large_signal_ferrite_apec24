function [s_mag, s_ele, s_tot] = get_post_power(E, J, H, dBdt)
% Compute the complex power.
%
%    Parameters:
%        E (float): electric field
%        J (float): current density
%        H (float): magnetic field
%        dBdt (float): magnetic flux density derivative
%
%    Returns:
%        s_mag (float): complex magnetic power
%        s_ele (float): complex electric power
%        s_tot (float): complex total power
%
%    Thomas Guillod - Dartmouth College.
%    2023 - MIT License.

s_mag = 0.5.*dBdt.*conj(H);
s_ele = 0.5.*E.*conj(J);
s_tot = s_mag+s_ele;

end