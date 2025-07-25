%% RC Battery Charging with Temperature & Non-Ideal Effects (DC)
clear; clc;

fprintf("========================================\n");
fprintf("  DC RC Battery Charging Simulation Tool\n");
fprintf("========================================\n");
fprintf("NOTE: This simulation models charging in a DC circuit using an RC network.\n");
fprintf("Units:\n");
fprintf("  - Capacitance: millifarads (mF)\n");
fprintf("  - Current limit: milliamps (mA)\n\n");

% User Input
getValidInput = @(prompt, minVal, maxVal) ...
    inputWithValidation(prompt, minVal, maxVal);

Vmax = getValidInput('Enter max battery voltage (V) [Recommended: 50–500, Max: 1000]: ', 0.1, 1000);
R0 = getValidInput('Enter nominal internal resistance at room temp (Ohms) [Recommended: 0.01–1.0, Max: 10]: ', 0.001, 10);
C_mF = getValidInput('Enter battery capacitance (millifarads) [Recommended: 1–5000, Max: 10000]: ', 0.1, 1e4);
I_limit_mA = getValidInput('Enter charging current limit (milliamps) [Recommended: 100–1000, Max: 10000]: ', 1, 1e4);
t_max = getValidInput('Enter total simulation time (seconds) [Recommended: 100–200, Max: 1000]: ', 100, 1000);
dt = getValidInput('Enter time step (e.g., 1): ', 0.0001, 1e4);
T = getValidInput('Enter ambient temperature (°C) [0–60]: ', -50, 150);

C = C_mF / 1000; 
charge_current_limit = I_limit_mA / 1000;

T0 = 25;               
alpha = 0.004;         
R = R0 * (1 + alpha * (T - T0));
fprintf('\nAdjusted Resistance at %.1f°C: %.4f Ohms\n\n', T, R);

t = 0:dt:t_max;
V = Vmax * (1 - exp(-t / (R * C)));
I = (Vmax - V) ./ R;
I = min(I, charge_current_limit);  % Apply current limit

P_input = Vmax * I;
E_input = trapz(t, P_input);  % Total energy supplied

P_loss_resistance = I.^2 * R;
E_loss_resistance = trapz(t, P_loss_resistance);

parasitic_loss_fraction = 0.05;  % 5%
E_loss_parasitic = parasitic_loss_fraction * E_input;

E_delivered = E_input - E_loss_resistance - E_loss_parasitic;

target_80 = 0.8 * Vmax;
target_100 = 0.999 * Vmax;

idx_80 = find(V >= target_80, 1);
idx_100 = find(V >= target_100, 1);

if isempty(idx_80)
    t_80 = NaN;
    warning('Could not determine time to reach 80%% charge.');
else
    t_80 = t(idx_80);
end

if isempty(idx_100)
    t_100 = NaN;
    warning('Could not determine time to reach ~100%% charge.');
else
    t_100 = t(idx_100);
end


%RESULTS 
fprintf("---------- Analysis Results ----------\n");
fprintf("Total Energy Supplied (from source): %.2f J\n", E_input);
fprintf("Energy Delivered to Battery (net): %.2f J\n", E_delivered);
fprintf("Energy Lost to Resistance: %.2f J\n", E_loss_resistance);
fprintf("Energy Lost to Parasitics: %.2f J\n", E_loss_parasitic);
fprintf("Time to 80%% charge: %.2f s\n", t_80);
fprintf("Time to ~100%% charge: %.2f s\n", t_100);
% PLOTS 
C_th = 100; 
T_batt = T + cumtrapz(t, P_loss_resistance) / C_th;
R_dynamic = R0 * (1 + alpha * (T_batt - T0));  

fprintf("\n========== Key Equations with Your Inputs ==========\n");

tau = R * C; 
fprintf("1. Voltage equation: V(t) = %.2f * (1 - exp(-t / %.4f))\n", Vmax, tau);
fprintf("2. Current equation: I(t) = (%.2f - V) / %.4f\n", Vmax, R);
fprintf("3. Power equation:   P(t) = %.2f * I(t)\n", Vmax);
fprintf("4. Temp-adjusted Resistance: R(T) = %.4f * (1 + %.4f * (%.1f - %.1f)) = %.4f Ohms\n", ...
    R0, alpha, T, T0, R);
fprintf("5. Battery Temp: T_batt(t) = %.1f + ∫(P_loss / %.1f) dt\n", T, C_th);
fprintf("RC Time Constant (τ): %.4f s\n", tau);
fprintf("====================================================\n\n");

figure;
plot(t, V, 'b', 'LineWidth', 2); grid on; hold on;
if ~isnan(t_80), xline(t_80, '--g', '80% Vmax'); end
if ~isnan(t_100), xline(t_100, '--r', '100% Vmax'); end
title('Voltage vs Time'); xlabel('Time [s]'); ylabel('Voltage [V]');
legend('Voltage', '80% Vmax', '100% Vmax');
grid on;

figure;
plot(t, I, 'r', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Current (A)');
title('Charging Current vs Time');
grid on;

figure;
plot(t, P_input, 'k', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Power (W)');
title('Charging Power vs Time');
grid on;

figure;
plot(t, T_batt, 'm', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Battery Temperature (°C)');
title('Battery Temperature vs Time');
grid on;

figure;
plot(t, R_dynamic, 'c', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Internal Resistance (Ohms)');
title('Internal Resistance Rise vs Time');
grid on;

temp_limit = 80; % max safe temp
idx_safe = find(T_batt < temp_limit, 1, 'last');
if ~isempty(idx_safe) && idx_safe < length(t)
    fprintf('⚠️  Simulation stopped early at t = %.2f s due to thermal limit (%.1f°C).\n', t(idx_safe), T_batt(idx_safe));
end


% INPUT VALIDATION 
function val = inputWithValidation(prompt, minVal, maxVal)
    while true
        val = input(prompt);
        if isempty(val) || ~isnumeric(val) || ~isscalar(val) || isnan(val)
            fprintf('❌ Please enter a numeric value.\n');
        elseif val <= minVal
            fprintf('❌ Value must be greater than %.2f.\n', minVal);
        elseif val > maxVal
            fprintf('⚠️  Value is above the recommended maximum of %.2f. Proceeding anyway...\n', maxVal);
            break;
        else
            break;
        end
    end
end
