%--------------------------------------------------------------------------
% This script simulates Fe-Mg diffusion profiles in clinopyroxene (Cpx)
% using a growth-then-diffusion approach. Ca profiles (assumed to reflect
% original growth zoning) are used to constrain the starting shape for 
% Fe-Mg diffusion.
%
% This code reproduces results comparable to Wang et al. (2023, Icarus),
% Figure 3. It estimates cooling rates and simulates profiles at
% different rates to visualize their impact on Fe-Mg diffusion.
% Author: Zilong Wang, written using MATLAB R2024b on 22th April, 2025
%--------------------------------------------------------------------------

clear; clc;

% -------------------------------------------------------------------------
% Load clinopyroxene profile data for five crystals
% -------------------------------------------------------------------------
cpx1 = readtable('CE5_cpx.xlsx','Sheet','cpx1','VariableNamingRule','preserve');
cpx2 = readtable('CE5_cpx.xlsx','Sheet','cpx2','VariableNamingRule','preserve');
cpx3 = readtable('CE5_cpx.xlsx','Sheet','cpx3','VariableNamingRule','preserve');
cpx4 = readtable('CE5_cpx.xlsx','Sheet','cpx4','VariableNamingRule','preserve');
cpx5 = readtable('CE5_cpx.xlsx','Sheet','cpx5','VariableNamingRule','preserve');

% Initial and final temperatures of the cooling path (°C)
T1 = 1033;
T2 = 950;

% -------------------------------------------------------------------------
% Fit Ca profiles to derive initial diffusion times for Fe-Mg
% -------------------------------------------------------------------------
[t_ini_1, ~, ~] = fit_to_Ca(cpx1, T1, T2, 0);
[t_ini_2, ~, ~] = fit_to_Ca(cpx2, T1, T2, 0);
[t_ini_3, ~, ~] = fit_to_Ca(cpx3, T1, T2, 0);
[t_ini_4, ~, ~] = fit_to_Ca(cpx4, T1, T2, 0);
[t_ini_5, ~, ~] = fit_to_Ca(cpx5, T1, T2, 0);

% -------------------------------------------------------------------------
% Cooling rates (K/hr) for bracketing the best-fit curve
% l = lower bound; m = best-fit; h = upper bound
% -------------------------------------------------------------------------
CR_l = [0.000015, 0.0001, 0.00005, 0.0002, 0.0001];
CR_m = [0.000024, 0.00024, 0.000066, 0.000302, 0.00011];
CR_h = [0.000035, 0.001, 0.0001, 0.001, 0.0002];
t_ini = [t_ini_1, t_ini_2, t_ini_3, t_ini_4, t_ini_5];

% -------------------------------------------------------------------------
% Generate profile comparison plots for each crystal and save high-res PNG
% -------------------------------------------------------------------------
for i = 1:5
    cpx = eval(['cpx', num2str(i)]);  % Access cpx1, cpx2, ...
    plot_diffusion(cpx, T1, T2, t_ini(i), CR_l(i), CR_m(i), CR_h(i));

    figureHandle = get(groot, 'CurrentFigure');
    print(figureHandle, sprintf('cpx%d.png', i), '-dpng', '-r900');
    pause(1)
end

close all;

function plot_diffusion(cpx, T1, T2, t_ini, CR_l, CR_m, CR_h)
%--------------------------------------------------------------------------
% plot_diffusion - Simulates and visualizes Fe-Mg diffusion in clinopyroxene
%
% This function:
%   (1) Fits Ca profiles using Fe-Mg diffusion coefficients to infer an
%       equivalent diffusion time (growth-then-diffusion initial condition).
%   (2) Simulates Fe-Mg diffusion at three cooling rates: lower, middle, and upper.
%   (3) Visualizes the results and compares simulated profiles with observed Mg#.
%
% INPUTS:
%   cpx     - Table containing compositional and spatial profile information
%   T1, T2  - Initial and final cooling temperatures (°C)
%   CR_ini  - Equivalent cooling rate derived from Ca profile (K/hr)
%   CR_l    - Lower bound of cooling rate for comparison (K/hr)
%   CR_m    - Best-estimate cooling rate (K/hr)
%   CR_h    - Upper bound of cooling rate for comparison (K/hr)
%
% OUTPUT:
%   Generates the 5 pannels of figure 3 in Wang et al. (2023)
%   Please see https://doi.org/10.1016/j.icarus.2022.115406 for the article
%--------------------------------------------------------------------------

% Convert cooling rate to diffusion time in years
t_l   = (T1 - T2) / (CR_l   * 24 * 365);
t_m   = (T1 - T2) / (CR_m   * 24 * 365);
t_h   = (T1 - T2) / (CR_h   * 24 * 365);

% Convert best-fit cooling rate (excluding Ca profile shaping) for labeling
cooling_rate_m = (T1 - T2) / ((t_m - t_ini) * 24 * 365); 

rate_str_l = format_rate_label(CR_l);
rate_str_m = format_rate_label(cooling_rate_m);
rate_str_h = format_rate_label(CR_h);

% Extract model geometry and boundary compositions for Ca and Fe-Mg
r_core_Ca = cpx.r_core(1);   r_rim_Ca = cpx.r_rim(1);
C0_core_Ca = cpx.C0_core(1); C0_rim_Ca = cpx.C0_rim(1);

r_core_FeMg = cpx.r_core_1(1); r_rim_FeMg = cpx.r_rim_1(1);
C0_core_FeMg = cpx.C0_core_1(1); C0_rim_FeMg = cpx.C0_rim_1(1);

% Spatial step for finite-difference
dr = 1;

% Experimental Mg# profile
d = cpx.("Distance from the core");
cpx_mgn = cpx.("Mg#");

% Simulate Ca profile as proxy for initial Fe-Mg
[R_ini, C_ini] = FeMg_Cpx(t_ini, T1, T2, r_core_Ca, r_rim_Ca, dr, C0_core_Ca, C0_rim_Ca);

% Simulate Fe-Mg profiles at three cooling rates
[R_l, C_l] = FeMg_Cpx(t_l, T1, T2, r_core_FeMg, r_rim_FeMg, dr, C0_core_FeMg, C0_rim_FeMg);
[R_m, C_m] = FeMg_Cpx(t_m, T1, T2, r_core_FeMg, r_rim_FeMg, dr, C0_core_FeMg, C0_rim_FeMg);
[R_h, C_h] = FeMg_Cpx(t_h, T1, T2, r_core_FeMg, r_rim_FeMg, dr, C0_core_FeMg, C0_rim_FeMg);

% Rescale Ca-derived profile to match the magnitude and center of Fe-Mg
C_ini_scaled = rescale_and_shift_to_match(C_ini, C_m);

% -------------------------------------------------------------------------
% Plotting
% -------------------------------------------------------------------------
figure; hold on;

% Simulated profiles and observed data
h1 = plot(R_ini, C_ini_scaled, '--k', 'LineWidth', 0.75, 'DisplayName', 'Initial profile');
h2 = plot(R_l, C_l, '-.', 'Color', [189 215 238]/255, 'LineWidth', 1.25, 'DisplayName', rate_str_l);
h3 = plot(R_h, C_h, '-.', 'Color', [157 195 230]/255, 'LineWidth', 1.25, 'DisplayName', rate_str_h);
h4 = plot(R_m, C_m, '-',  'Color', [46 117 182]/255,  'LineWidth', 1.5,  'DisplayName', rate_str_m);
h0 = scatter(d, cpx_mgn, 40, 'k', 'Marker', 'd', 'DisplayName', 'Observed');

% Axis and annotation setup
xlabel('Distance from core (\mum)', 'FontName','Calibri', 'FontSize', 12);
ylabel('Mg#', 'FontName','Calibri', 'FontSize', 12);
legend([h3,h4,h2,h1,h0], 'Location', 'best', 'FontSize', 12, 'FontName', 'Calibri');

set(gca, 'Box', 'on', ...
         'FontName', 'Calibri', 'FontSize', 12, ...
         'LineWidth', 0.75, ...
         'TickDir', 'out', ...
         'TickLength', [0.01 0.01]);

end

%--------------------------------------------------------------------------
function rate_str = format_rate_label(rate)
% Formats cooling rate in K/h with adaptive precision for display
    if rate < 1e-4
        fmt = '%.6f';
    elseif rate < 1e-3
        fmt = '%.5f';
    elseif rate < 1e-2
        fmt = '%.4f';
    else
        fmt = '%.3f';
    end
    rate_str = sprintf([fmt ' K/h'], rate);
    rate_str = regexprep(rate_str, '0+$', '');  % remove trailing zeros
    rate_str = regexprep(rate_str, '\.$', '');  % remove trailing decimal
end

function C_ini_shifted = rescale_and_shift_to_match(C_ini, C_target)
% RESCALE_AND_SHIFT_TO_MATCH - Linearly rescales C_ini to match C_target's range,
% then shifts it so their centers align (based on max gradient), with edge filling.
%
% INPUT:
%   C_ini     - Original input array to transform (e.g., erf-like array)
%   C_target  - Reference array to match (for range and center)
%
% OUTPUT:
%   C_ini_shifted - Transformed array, rescaled and center-aligned

% 1. Rescale C_ini to match C_target range
minT = min(C_target);
maxT = max(C_target);

C_ini_scaled = rescale(C_ini, minT, maxT);  % Linearly scale to same range

% 2. Estimate center using max gradient (symmetric point)
[~, center_ini] = max(abs(gradient(C_ini_scaled)));
[~, center_target] = max(abs(gradient(C_target)));

% 3. Shift C_ini_scaled to align centers
shift = center_target - center_ini;
N = length(C_ini_scaled);
C_ini_shifted = zeros(size(C_ini_scaled));

if shift > 0
    % Shift right: pad head with left value
    C_ini_shifted(1:shift) = C_ini_scaled(1);
    C_ini_shifted(shift+1:end) = C_ini_scaled(1:end-shift);
elseif shift < 0
    % Shift left: pad tail with right value
    shift = abs(shift);
    C_ini_shifted(end-shift+1:end) = C_ini_scaled(end);
    C_ini_shifted(1:end-shift) = C_ini_scaled(shift+1:end);
else
    % No shift
    C_ini_shifted = C_ini_scaled;
end

end