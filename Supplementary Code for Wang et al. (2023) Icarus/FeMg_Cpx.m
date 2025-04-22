function [R, C] = FeMg_Cpx(year, T1, T2, r_core, r_rim, dr, C0_core, C0_rim)
% FEMG_CPX calculates Fe-Mg diffusion in clinopyroxene using finite difference.
% 
% This function simulates the 1D radial diffusion profile of Fe-Mg between core and rim 
% in a clinopyroxene grain over geological timescales using a temperature-dependent 
% diffusion coefficient and time-varying temperature.
% 
% INPUTS:
%   year      - total duration of diffusion (years)
%   T1        - initial temperature (Celsius)
%   T2        - final temperature (Celsius)
%   r_core    - core radius (microns or other length unit)
%   r_rim     - rim radius (microns or other length unit)
%   dr        - spatial step (same unit as r_core/r_rim)
%   C0_core   - initial Fe-Mg content in core
%   C0_rim    - initial Fe-Mg content in rim
% 
% OUTPUTS:
%   R         - radial positions (concatenated core and rim)
%   C         - Fe-Mg concentrations after diffusion

% Convert total diffusion time from years to seconds
t_total = year * 365 * 24 * 3600;

% Time step and number of steps
dt = year * 1e4;               % Fixed time increment (in seconds)
t_nstep = round(t_total / dt); % Total number of time steps

% Temperature change rate (C/s)
dT = (T1 - T2) / t_total;

% Spatial grids for core and rim
R_core = 0:dr:r_core;
R_rim  = r_core:dr:r_rim;

n_core = length(R_core);
n_rim  = length(R_rim);
n_total = n_core + n_rim - 1; % Total number of radial nodes

% Initialize concentrations
C_core = ones(n_core, 1) * C0_core;
C_rim  = ones(n_rim, 1)  * C0_rim;

% Preallocate arrays
B = zeros(n_total); % Tridiagonal matrix for implicit scheme
a = zeros(n_total, 1);
b = zeros(n_total, 1);
c = zeros(n_total, 1);
d = zeros(n_total, 1);
Cs = zeros(n_total, 1);

Tn = T1; % Initial temperature

% Time stepping loop
for k = 1:t_nstep
    % Update temperature linearly
    Tn = Tn - dT * dt;

    % Temperature-dependent diffusivity (Arrhenius law, D in m^2/s)
    D = 2.77e5 * exp(-320.7 / (0.008314 * (Tn + 273.15)));

    % Alpha = D*dt / (2*dr^2) for Crank-Nicolson scheme
    alpha = D * dt / (2 * dr^2);

    % Assemble coefficient matrix B and RHS vector d
    for i = 2:n_total-1
        a(i) = -alpha * (1 - dr / ((i-1) * dr));
        b(i) = 1 + 2 * alpha;
        c(i) = -alpha * (1 + dr / ((i-1) * dr));

        B(i, i-1) = -a(i);
        B(i, i)   = 1 - 2 * alpha;
        B(i, i+1) = -c(i);
    end

    % Boundary conditions
    % Inner boundary (Neumann)
    B(1, 1) = 1;
    B(1, 2) = -1;
    b(1) = -1;
    c(1) = 1;

    % Core-rim boundary (continuity)
    B(n_core, n_core-1) = -D;
    B(n_core, n_core)   = 2*D;
    B(n_core, n_core+1) = -D;
    a(n_core) = D;
    b(n_core) = -2*D;
    c(n_core) = D;

    % Outer boundary (Neumann)
    B(n_total, n_total-1) = 1;
    B(n_total, n_total)   = -1;
    a(n_total) = -1;
    b(n_total) = 1;

    % Construct current profile
    C_full = [C_core(1:n_core-1); C_rim];

    % Compute RHS
    d = B * C_full;

    % Solve tridiagonal system using Thomas algorithm
    for i = 2:n_total
        b(i) = b(i) - a(i) * c(i-1) / b(i-1);
        d(i) = d(i) - a(i) * d(i-1) / b(i-1);
    end
    Cs(n_total) = d(n_total) / b(n_total);
    for i = n_total-1:-1:1
        Cs(i) = (d(i) - c(i) * Cs(i+1)) / b(i);
    end

    % Update concentration profile
    C_core(1:n_core-1) = Cs(1:n_core-1);
    C_rim = Cs(n_core:end);

    % Ensure continuity at core-rim boundary
    C_avg = (C_core(n_core) + C_rim(1)) / 2;
    C_core(n_core) = C_avg;
    C_rim(1) = C_avg;
end

% Concatenate output radius and concentration
R = [R_core, R_rim(2:end)];
C = [C_core; C_rim(2:end)];

end