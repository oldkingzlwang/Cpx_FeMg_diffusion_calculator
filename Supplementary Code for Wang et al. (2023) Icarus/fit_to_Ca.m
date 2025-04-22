function [best_year, best_CR, min_misfit] = fit_to_Ca(cpx, T1, T2, plot_flag)
% FIT_TO_CA - Estimate diffusion time by fitting a synthetic Fe-Mg diffusion
%             profile to an observed Ca profile in clinopyroxene.
%
% This function assumes that the shape of Ca zoning (often unaffected by 
% diffusional modification) reflects the original Fe-Mg profile at the 
% onset of diffusion. It searches for the best-fit diffusion timescale 
% (in years) that produces a synthetic Fe-Mg profile morphologically 
% similar to the observed Ca profile.
%
% INPUT:
%   cpx       - Struct or table containing clinopyroxene profile data
%               with fields: Distance from core, Ca content, r_core, r_rim,
%               and initial Ca values at core and rim.
%   T1, T2    - Initial and final temperatures (in °C) for the cooling model
%   plot_flag - Boolean (true/false) to display comparison plots
%
% OUTPUT:
%   best_year - Best-fit diffusion time (years)
%   best_CR   - Estimated cooling rate (K/hr) over best_year
%   min_misfit - Minimum misfit value (variance between profiles)

% Extract distance and Ca zoning profile
d = cpx.("Distance from the core");
cpx_Ca = cpx.Ca;

% Set up geometry and initial conditions
dr = 1;  % Spatial resolution in microns
r_core = cpx.r_core(1);
r_rim  = cpx.r_rim(1);
C0_core = cpx.C0_core(1);
C0_rim  = cpx.C0_rim(1);

% Define search range for diffusion time (in years)
year_range = linspace(0.01, 100, 100);  % Log-linear or linear range

% Initialize misfit array
misfit = zeros(size(year_range));

% Loop over candidate diffusion durations
for i = 1:length(year_range)
    year = year_range(i);

    % Generate Fe-Mg diffusion profile under given cooling scenario
    [R, C] = FeMg_Cpx(year, T1, T2, r_core, r_rim, dr, C0_core, C0_rim);

    % Interpolate modeled profile to observed positions
    C_interp = interp1(R, C, d, 'linear', 'extrap');

    % Calculate misfit as variance between modeled and observed profiles
    misfit(i) = var(C_interp - cpx_Ca);
end

% Extract best-fit diffusion time and compute corresponding cooling rate
[min_misfit, idx] = min(misfit);
best_year = year_range(idx);
best_CR = abs(T1 - T2) / (best_year * 365 * 24);  % K/hr

% Plot results if requested
if plot_flag
    % Plot profile fit
    figure;
    [R, C] = FeMg_Cpx(best_year, T1, T2, r_core, r_rim, dr, C0_core, C0_rim);
    hold on;
    scatter(d, cpx_Ca, 30, 'filled', 'd', 'DisplayName', 'Observed Ca');
    plot(R, C, 'LineWidth', 1.5, 'DisplayName', 'Best-fit Fe-Mg');
    hold off;
    xlabel('Distance from core (μm)');
    ylabel('Ca content (apfu)');
    legend('Location', 'best');
    title(['Best-fit diffusion time: ', num2str(best_year, '%.2f'), ' yrs']);

    % Plot misfit evolution over time
    figure;
    plot(year_range, misfit, 'k-', 'LineWidth', 1.5);
    xlabel('Diffusion time (years)');
    ylabel('Misfit (variance)');
    title('Misfit between modeled and observed profiles');
end

end