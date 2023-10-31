function run_single()
% Solve the magnetic-dielectric effects in a core with a given cross-section.
%
%    The quasi-static approximation and the full wave solution are compared.
%    The large-signal material parameters are used for the permeability and permitivitty.
%
%    Thomas Guillod - Dartmouth College.
%    2023 - MIT License.

close('all')
addpath(genpath('code'))

% define the parameters
f = 300e3;
A = 10e-4;
B_src = 100e-3;

% get the solver parameters
[problem, param, mat] = get_param(f, B_src, A);

% solve the quasi-static approximation
[sol_approx, is_valid] = get_solve_approx(problem, mat);
assert(is_valid, 'extrapolation for approx solution')

% solve the full wave solution
[sol_wave, is_valid] = get_solve_wave(problem, param, mat);
assert(is_valid, 'extrapolation for wave solution')

% display the scalar parameters
fprintf('E\n')
fprintf('    approx = %.3f V/mm\n', 1e-3.*sol_approx.E)
fprintf('    wave = %.3f V/mm\n', 1e-3.*sol_wave.E)
fprintf('D\n')
fprintf('    approx = %.3f nC/mm2\n', 1e3.*sol_approx.D)
fprintf('    wave = %.3f nC/mm2\n', 1e3.*sol_wave.D)
fprintf('J\n')
fprintf('    approx = %.3f mA/mm2\n', 1e-3.*sol_approx.J)
fprintf('    wave = %.3f mA/mm2\n', 1e-3.*sol_wave.J)
fprintf('B\n')
fprintf('    approx = %.3f mT\n', 1e3.*sol_approx.B)
fprintf('    wave = %.3f mT\n', 1e3.*sol_wave.B)
fprintf('H\n')
fprintf('    approx = %.3f mA/mm\n', 1e0.*sol_approx.H)
fprintf('    wave = %.3f mA/mm\n', 1e0.*sol_wave.H)
fprintf('s_mag\n')
fprintf('    approx = %.2f + %.2fj mVA/cm3\n', 1e-3.*real(sol_approx.s_mag), +1e-3.*imag(sol_approx.s_mag))
fprintf('    wave = %.2f + %.2fj mVA/cm3\n', 1e-3.*real(sol_wave.s_mag), +1e-3.*imag(sol_wave.s_mag))
fprintf('s_ele\n')
fprintf('    approx = %.2f - %.2fj mVA/cm3\n', 1e-3.*real(sol_approx.s_ele), -1e-3.*imag(sol_approx.s_ele))
fprintf('    wave = %.2f - %.2fj mVA/cm3\n', 1e-3.*real(sol_wave.s_ele), -1e-3.*imag(sol_wave.s_ele))
fprintf('s_tot\n')
fprintf('    approx = %.2f + %.2fj mVA/cm3\n', 1e-3.*real(sol_approx.s_tot), +1e-3.*imag(sol_approx.s_tot))
fprintf('    wave = %.2f + %.2fj mVA/cm3\n', 1e-3.*real(sol_wave.s_tot), +1e-3.*imag(sol_wave.s_tot))
fprintf('mu\n')
fprintf('    approx = %.2f - %.2fj p.u.\n', 1e0.*real(sol_approx.mu), -1e0.*imag(sol_approx.mu))
fprintf('    wave = %.2f - %.2fj p.u.\n', 1e0.*real(sol_wave.mu), -1e0.*imag(sol_wave.mu))
fprintf('eps\n')
fprintf('    approx = %.2f - %.2fj p.u.\n', 1e0.*real(sol_approx.eps), -1e0.*imag(sol_approx.eps))
fprintf('    wave = %.2f - %.2fj p.u.\n', 1e0.*real(sol_wave.eps), -1e0.*imag(sol_wave.eps))
fprintf('lambda\n')
fprintf('    approx = %.3f mm\n', 1e3.*sol_approx.lambda)
fprintf('    wave = %.3f mm\n', 1e3.*sol_wave.lambda)
fprintf('delta\n')
fprintf('    approx = %.3f mm\n', 1e3.*sol_approx.delta)
fprintf('    wave = %.3f mm\n', 1e3.*sol_wave.delta)

% plot the magnetic flux density and the electric field
figure()
subplot(2,1,1)
plot(1e3.*sol_wave.r_vec, 1e3.*abs(sol_wave.B_vec), 'LineWidth', 2.0)
grid('on')
xlabel('r (mm)')
ylabel('B (mT)')
title('Mangetic Flux Density')
subplot(2,1,2)
plot(1e3.*sol_wave.r_vec, 1e-3.*abs(sol_wave.E_vec), 'LineWidth', 2.0)
grid('on')
xlabel('r (mm)')
ylabel('E (V/mm)')
title('Electric Field')

end
