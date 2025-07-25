%% BMW i7 DC Fast‑Charge Simulation — MATLAB version

clear; clc;

E_kwh   = 101.7;                 
E_target = E_kwh * 3.6e6;        

P_const = 150e3;                 
t_total = E_target / P_const;     

t  = linspace(0, t_total, 2001)'; 
dt = t(2) - t(1);

V0 = 290;  V_nom = 376;
tau = t_total / 3;                
V = V_nom - (V_nom - V0).*exp(-t/tau);

I = P_const ./ V;

P_curve = P_const * ones(size(t));

% Numerics
E_delivered = trapz(t, P_curve);            

R_int = 0.03;                               
Q_loss = trapz(t, I.^2 * R_int);           

dVdt = gradient(V, dt);                     

E_cum   = cumsum(P_curve) * dt;             
idx_80  = find(E_cum >= 0.80*E_target, 1);
time_80 = t(idx_80)/60;                     
time_100 = t_total/60;                      
time_80_100 = time_100 - time_80;           

seg_len  = floor(numel(t)/3);
avg_rate = zeros(1,3);
for k = 1:3
    if k < 3
        seg_idx = (k-1)*seg_len + (1:seg_len);
    else
        seg_idx = (k-1)*seg_len + 1 : numel(t);
    end
    avg_rate(k) = mean(dVdt(seg_idx))*60;   % V/min
end

E_useful = E_delivered - Q_loss;

% Console output 
fprintf('=== BMW i7 DC Fast‑Charge Simulation Results ===\n');
fprintf('State of Charge = SOC\n');
fprintf('Total energy delivered (input)    : %6.2f kWh (%.2f MJ)\n', ...
        E_delivered/3.6e6, E_delivered/1e6);
fprintf('Energy lost to Resistance (loss)    : %6.2f kWh (%.2f MJ, %4.1f %% of input)\n', ...
        Q_loss/3.6e6, Q_loss/1e6, 100*Q_loss/E_delivered);
fprintf('Total useful energy delivered (output): %6.2f kWh (%.2f MJ)\n', ...
        E_useful/3.6e6, E_useful/1e6);
fprintf('Charging time to 80 %%               : %6.2f min\n', time_80);
fprintf('Charging time to 100 %%              : %6.2f min\n', time_100);
fprintf('Time 80 %% → 100 %% SOC              : %6.2f min\n', time_80_100);

fprintf('Average dV/dt (0‑33%%)               : %6.2f V/min\n', avg_rate(1)); 
fprintf('Average dV/dt (33‑67%%)              : %6.2f V/min\n', avg_rate(2)); 
fprintf('Average dV/dt (67‑100%%)             : %6.2f V/min\n', avg_rate(3));
overallAvgRate = mean(dVdt) * 60;  % V/min
fprintf('Average dV/dt (overall)             : %6.2f V/min\n\n', overallAvgRate);

time_min = t/60;
xLim = [0, time_min(end)];

% Plot

figure(1);
plot(time_min, V, 'LineWidth',1.5); hold on;
plot(time_80,  V(idx_80), 'ro','MarkerFaceColor','r');
plot(time_100, V(end)  , 'go','MarkerFaceColor','g');
text(time_80  +1, V(idx_80), '80%','Color','r');
text(time_100 +1, V(end)  , '100%','Color','g');
xlabel('Time [min]'); ylabel('Voltage [V]');
title('Voltage vs. Time');
xlim(xLim);
grid on; ax = gca; ax.XMinorGrid='on'; ax.YMinorGrid='on';
legend('Voltage','80% SOC','100% SOC','Location','best'); hold off;

figure(2);
plot(time_min, I, 'LineWidth',1.5); hold on;
plot(time_80 , I(idx_80), 'ro','MarkerFaceColor','r');
plot(time_100, I(end)  , 'go','MarkerFaceColor','g');
text(time_80  +1, I(idx_80), '80%','Color','r');
text(time_100 +1, I(end)  , '100%','Color','g');
xlabel('Time [min]'); ylabel('Current [A]');
title('Current vs. Time');
xlim(xLim);
grid on; ax = gca; ax.XMinorGrid='on'; ax.YMinorGrid='on';
legend('Current','80% SOC','100% SOC','Location','best'); hold off;

figure(3);
plot(time_min, P_curve/1000, 'LineWidth',1.5); hold on;
xlabel('Time [min]'); ylabel('Power [kW]');
title('Power vs. Time');
xlim(xLim);
grid on; ax = gca; ax.XMinorGrid='on'; ax.YMinorGrid='on';
yl = ylim;
patch([time_80 time_100 time_100 time_80], ...
      [yl(1) yl(1) yl(2) yl(2)], ...
      [0.85 0.85 0.92], 'FaceAlpha',0.25, 'EdgeColor','none');
xline(time_80 , '--', ' 80%','LabelVerticalAlignment','bottom');
xline(time_100, '--', '100%','LabelVerticalAlignment','bottom');

text(time_100 - 0.1*diff(xLim), yl(2)*0.95, ...
     sprintf('Total Energy = %.2f kWh', E_delivered/3.6e6), ...
     'FontSize',9, 'BackgroundColor','white');

legend('Power','80–100 % window','Location','best'); hold off;

figure(4);
plot(time_min, dVdt*60, 'b', 'LineWidth',1.5);   
xlabel('Time [min]'); ylabel('dV/dt [V/min]');
title('Rate of Voltage Change (dV/dt) vs. Time');
xlim(xLim);
grid on; ax = gca; ax.XMinorGrid='on'; ax.YMinorGrid='on';

% Advanced Temperature
T_amb     = 25;          
T_init    = 25;          
C_th      = 180e3;       
R_th      = 0.35;        
alpha_R   = 0.04;        
P_min     = 50e3;        

n  = numel(t);
T  = zeros(n,1);    T(1)  = T_init;      
R  = zeros(n,1);    R(1)  = R_int;        
P  = zeros(n,1);                      
I2 = zeros(n,1);                      
E_cum = 0;                            
soc   = zeros(n,1);

for k = 1:n
    soc(k) = E_cum / E_target;
    if soc(k) < 0.80
        P(k) = P_const;
    else
        frac = (soc(k) - 0.80) / 0.20;           % 0 → 1
        P(k) = P_const - frac*(P_const - P_min);
    end

    I2(k) = P(k) / V(k);         
    P_loss_k = I2(k)^2 * R(k);
    if k < n
        dTdt = (P_loss_k - (T(k)-T_amb)/R_th) / C_th;
        T(k+1) = T(k) + dTdt*dt;
    end
    if k < n
        E_cum = E_cum + P(k)*dt;
    end
    if k < n
        R(k+1) = R_int * exp(alpha_R*(25 - T(k+1)));
    end
    if soc(k) >= 1.0
        n = k;            % truncate arrays neatly
        t   = t(1:n);     time_min = t/60;
        V   = V(1:n);     R   = R(1:n);
        T   = T(1:n);     P   = P(1:n);
        I2  = I2(1:n);    soc = soc(1:n);
        break
    end
end

time_full_T = t(end)/60;                  
E_deliv_T   = trapz(t, P);                   
Q_loss_T    = trapz(t, I2.^2 .* R);         
E_useful_T  = E_deliv_T - Q_loss_T;

fprintf('\n=== Temperature & Non‑Ideal Model Results ================\n');
fprintf('Total charging time (0–100%%)       : %6.2f min\n', time_full_T);
fprintf('Total energy delivered (input)     : %6.2f kWh (%.2f MJ)\n', ...
        E_deliv_T/3.6e6, E_deliv_T/1e6);
fprintf('Energy lost to R_int(T)            : %6.2f kWh (%.2f MJ, %4.1f %% of input)\n', ...
        Q_loss_T/3.6e6, Q_loss_T/1e6, 100*Q_loss_T/E_deliv_T);
fprintf('Pack temperature rise (ΔT)          : %6.2f °C (%.1f → %.1f °C)\n', ...
        T(end)-T_init, T_init, T(end));
fprintf('-----------------------------------------------------------\n');

% PLOTS 
figure(5);
plot(time_min, T, 'r','LineWidth',1.6); grid on;
xlabel('Time [min]'); ylabel('Temperature [°C]');
title('Battery Temperature vs. Time');

figure(6);
plot(time_min, R*1e3, 'k','LineWidth',1.6); grid on;
xlabel('Time [min]'); ylabel('R_{int}(T) [mΩ]');
title('Internal Resistance Rise vs. Time');
