function [sol, is_valid, iter] = get_solve_wave(problem, param, mat)
% Solve the magnetic-dielectric effects in a core with a full wave solution.
%
%    A infinite cylindrical core is considered.
%    A full wave solution is used (with wave propagation).
%    The (average) flux density through the cross-section is imposed (excitation).
%    The material parameters can be amplitude-dependent (locally linearized).
%    The material parameters are variable over the cross-section (spatially dependent).
%
%    The magnetic field is flowing in the out-of-plane direction (z-axis).
%    The magnetic field is only dependent on the radial position (rho-axis).
%
%    The electric field is flowing in the in-plane direction (phi-axis).
%    The electric field is only dependent on the radial position (rho-axis).
%
%    Parameters:
%        problem (struct): description of the problem (frequency, flux density, radius)
%        param (struct): description of the solver numerical parameters and tolerances
%        mat (struct): description of the material parameters (permeability, permitivitty)
%
%    Returns:
%        sol (struct): solution of the full wave electromagnetic problem
%        is_valid (boolean): validity of the solution (detect extrapolation)
%        iter (struct): information about the iterative solver convergence
%
%    Thomas Guillod - Dartmouth College.
%    2023 - MIT License.

% extract solver parameters
iter_min = param.iter_min;
iter_tol = param.iter_tol;
iter_relax = param.iter_relax;

% get the spatially dependent material parameters
%     - for the first iteration, the material parameters are constant
%     - for the first iteration, the quasi-static approximation is used
space = get_material_init(problem, mat);

% solve the electromagnetic problem
[sol, is_valid] = get_solve_single(problem, param, mat, space);

% init the iteration
iter_nb = 0;
iter_stop = false;

% iterative solver
%    - update the material parameters
%    - solve the electromagnetic problem
%    - check for convergence
while iter_stop==false
    % get the spatially dependent material parameters
    %     - the material parameters are spatially dependent
    %     - the material parameters are updated with a relaxation parameter
    space = get_material_iter(sol, mat, space, iter_relax);

    % solve the electromagnetic problem
    [sol_new, is_valid] = get_solve_single(problem, param, mat, space);

    % check the error (on the complex power) between the last two iterations
    [err_mag, err_ele, iter_success] = get_convergence(sol, sol_new, iter_tol);

    % update the iteration counter
    iter_nb = iter_nb+1;

    % update the solution for the next iteration
    sol = sol_new;

    % check if convergence is achieved
    iter_stop = (iter_nb>=iter_min)&(iter_success==true);

    % information about the iterative solver
    iter = struct('iter_nb', iter_nb, 'err_ele', err_ele, 'err_mag', err_mag);
end

end
