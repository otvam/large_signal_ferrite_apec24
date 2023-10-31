function [err_mag, err_ele, iter_success] = get_convergence(sol, sol_new, iter_tol)
% Check the convergence (on the complex power) between two iterations.
%
%    Parameters:
%        sol (struct): existing solution of the full wave electromagnetic problem
%        sol_new (struct): new solution of the full wave electromagnetic problem
%        iter_tol (float): relative tolerance for convergence
%
%    Returns:
%        err_mag (float): relative error on the magnetic power
%        err_ele (float): relative error on the electric power
%        iter_success (boolean): flag indicating if the convergence is achieved
%
%    Thomas Guillod - Dartmouth College.
%    2023 - MIT License.

% extract the solution
s_mag = sol.s_mag;
s_ele = sol.s_ele;
s_mag_new = sol_new.s_mag;
s_ele_new = sol_new.s_ele;

% compute relative error
err_mag = abs(s_mag-s_mag_new)./abs(s_mag_new);
err_ele = abs(s_ele-s_ele_new)./abs(s_ele_new);

% check for convergence
iter_success = (err_mag<=iter_tol)&(err_ele<=iter_tol);

end
