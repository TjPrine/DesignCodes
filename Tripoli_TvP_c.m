%% Setup
clc;
clear;

%% Initial Conditions
Tank_total_Pressure_init = 3.1e+7; %Pa
Tank_total_Temp_init = 298.15; %K
Thrust_eff = 0.90; % Efficiency of Nozzle
P_atm = 1.01e5; %Pa - Exit Pressure Equals Atm pressure
d_star = 4.0e-3; %m
%d_exit = 2.5e-2; %m
P_c = .01*Tank_total_Pressure_init; %Pa
k = 1.4; % Specific Heat Ratio
R_gas = 296; %J/kg/K\

%m_tank_init = Tank_total_Pressure_init*V_tank/R_gas/Tank_total_Temp_init;
A_star = pi/4*d_star^2;

%% Simulate Change in Thrust as Chamber Pressure Increases
P_chamber(1) = P_c;
n = 1;
pressure_step = .01*Tank_total_Pressure_init;
while P_chamber(n) < Tank_total_Pressure_init

    C_f(n) = sqrt(2*k*k/(k-1)*(2/(k+1))^((k+1)/(k-1))*(1 - (P_atm/P_chamber(n))^((k-1)/k)));
    Thrust(n) = P_chamber(n)*A_star*C_f(n)*Thrust_eff;
    
    %Update iteration values
    P_chamber(n+1) = P_chamber(n) + pressure_step;
    n = n+1; 
end
n = n-1;
figure(1)
plot(P_chamber(1:n),Thrust(1:n))
title('Thrust v Chamber Pressure')
xlabel('Chamber Pressure(Pa)');
ylabel('Thrust(N)');
