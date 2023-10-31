function run_sweep()
% Solve the magnetic-dielectric effects in cores with different cross-sections.
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
B_src = 100e-3;
A = logspace(log10(0.1e-4), log10(10.0e-4), 30);

% solve the electromagnetic problems
for i=1:length(A)
    % show solver progress
    fprintf('%d / %d\n', i, length(A))

    % get the solver parameters
    [problem, param, mat, s_ref_tmp] = get_param(f, B_src, A(i));

    % assign the complex power without dielectric effects
    s_ref(i) = s_ref_tmp;

    % solve the quasi-static approximation
    [sol_approx{i}, is_valid] = get_solve_approx(problem, mat);
    assert(is_valid, 'extrapolation for approx solution')

    % solve the full wave solution
    [sol_wave{i}, is_valid] = get_solve_wave(problem, param, mat);
    assert(is_valid, 'extrapolation for wave solution')
end

% extract the scalar variables (complex powers and fields)
for i=1:length(A)
    s_mag_approx(i) = sol_approx{i}.s_mag;
    s_ele_approx(i) = sol_approx{i}.s_ele;
    s_tot_approx(i) = sol_approx{i}.s_tot;
    B_approx(i) = sol_approx{i}.B;
    E_approx(i) = sol_approx{i}.E;

    s_mag_wave(i) = sol_wave{i}.s_mag;
    s_ele_wave(i) = sol_wave{i}.s_ele;
    s_tot_wave(i) = sol_wave{i}.s_tot;
    B_wave(i) = sol_wave{i}.B;
    E_wave(i) = sol_wave{i}.E;
end

% compute the impact of the dielectric effects on the complex power
fact_P_approx = real(s_tot_approx)./real(s_ref);
fact_Q_approx = imag(s_tot_approx)./imag(s_ref);
fact_P_wave = real(s_tot_wave)./real(s_ref);
fact_Q_wave = imag(s_tot_wave)./imag(s_ref);

% plot the magnetic flux density and the electric field RMS values
figure()
subplot(2,1,1)
semilogx(1e4.*A, 1e3.*B_wave, '-', 'LineWidth', 2.0)
hold('on')
semilogx(1e4.*A, 1e3.*B_approx, '--', 'LineWidth', 2.0)
grid('on')
legend('wave', 'approx')
xlabel('A (cm2)')
ylabel('B (mT)')
title('Magnetic Flux Density')
subplot(2,1,2)
semilogx(1e4.*A, 1e-3.*E_wave, '-', 'LineWidth', 2.0)
hold('on')
semilogx(1e4.*A, 1e-3.*E_approx, '--', 'LineWidth', 2.0)
grid('on')
legend('wave', 'approx')
xlabel('A (cm2)')
ylabel('E (E/mm)')
title('Magnetic Field')

% plot the complex power average values
figure()
subplot(2,1,1)
semilogx(1e4.*A, 1e-3.*real(s_mag_wave), '-', 'LineWidth', 2.0)
hold('on')
semilogx(1e4.*A, 1e-3.*real(s_mag_approx), '--', 'LineWidth', 2.0)
semilogx(1e4.*A, 1e-3.*real(s_ele_wave), '-', 'LineWidth', 2.0)
semilogx(1e4.*A, 1e-3.*real(s_ele_approx), '--', 'LineWidth', 2.0)
semilogx(1e4.*A, 1e-3.*real(s_tot_wave), '-', 'LineWidth', 2.0)
semilogx(1e4.*A, 1e-3.*real(s_tot_approx), '--', 'LineWidth', 2.0)
grid('on')
legend('mag / wave', 'mag / approx', 'ele / wave', 'ele / approx', 'tot / wave', 'tot / approx')
xlabel('A (cm2)')
ylabel('P (mW/cm3)')
title('Active Power')
subplot(2,1,2)
semilogx(1e4.*A, +1e-3.*imag(s_mag_wave), '-', 'LineWidth', 2.0)
hold('on')
semilogx(1e4.*A, +1e-3.*imag(s_mag_approx), '--', 'LineWidth', 2.0)
semilogx(1e4.*A, -1e-3.*imag(s_ele_wave), '-', 'LineWidth', 2.0)
semilogx(1e4.*A, -1e-3.*imag(s_ele_approx), '--', 'LineWidth', 2.0)
semilogx(1e4.*A, +1e-3.*imag(s_tot_wave), '-', 'LineWidth', 2.0)
semilogx(1e4.*A, +1e-3.*imag(s_tot_approx), '--', 'LineWidth', 2.0)
grid('on')
legend('mag / wave', 'mag / approx', 'ele / wave', 'ele / approx', 'tot / wave', 'tot / approx')
xlabel('A (cm2)')
ylabel('P (mVA/cm3)')
title('Reactive Power')

% plot the impact of the dielectric effects on the complex power
figure()
subplot(2,1,1)
semilogx(1e4.*A, fact_P_approx, '--', 'LineWidth', 2.0)
hold('on')
semilogx(1e4.*A, fact_P_wave, '-', 'LineWidth', 2.0)
grid('on')
legend('approx', 'wave')
xlabel('A (cm2)')
ylabel('factor (p.u.)')
title('Ratio / Active Power')
subplot(2,1,2)
semilogx(1e4.*A, fact_Q_approx, '--', 'LineWidth', 2.0)
hold('on')
semilogx(1e4.*A, fact_Q_wave, '-', 'LineWidth', 2.0)
grid('on')
legend('approx', 'wave')
xlabel('A (cm2)')
ylabel('factor (p.u.)')
title('Ratio / Reactive Power')

end
