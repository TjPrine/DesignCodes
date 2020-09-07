% About:
%   Cold Gas Thruster Project
%   Author: Sean Tibbetts
%   Purpose: Simulate the flow of a Cold-Gas Thruster

%% Setup
clc;
clear;

%% Initial Conditions
Tank_total_Pressure_init = 2.9303e+7; %Pa
Tank_total_Temp_init = 298.15; %K
Thrust_eff = 0.98; % Efficiency of Nozzle
P_atm = 1e5; %Pa - Exit Pressure Equals Atm pressure
d_star = 4.0e-3; %m
%d_exit = 2.5e-2; %m
P_c = 2.06e6; %Pa
V_tank = 0.00147484/2; %m^3
k = 1.4; % Specific Heat Ratio
R_gas = 287; %J/kg/K\

m_tank_init = Tank_total_Pressure_init*V_tank/R_gas/Tank_total_Temp_init;
A_star = pi/4*d_star^2;

%% Simulate the CGT using Total Temp and Total Pres Change
m_expelled(1) = 0.0;
P_tank(1) = Tank_total_Pressure_init;
T_tank(1) = Tank_total_Temp_init;
P_chamber(1) = P_c;
m_tank(1) = m_tank_init;
n = 1;
t_step = 0.1;
time = 0;
while P_tank(n) > 101325
    T_chamber(n) = T_tank(n);%*(P_chamber(n)/P_tank(n))^((k-1)/k);
    mdot(n) = P_chamber(n)*A_star*sqrt(k/R_gas/T_chamber(n))*(2/(k+1))^((k+1)/2/(k-1));
    m_expelled(n+1) = m_expelled(n) + t_step*mdot(n);
    C_f(n) = sqrt(2*k*k/(k-1)*(2/(k+1))^((k+1)/(k-1))*(1 - (P_atm/P_chamber(n))^((k-1)/k)));
    Thrust(n) = P_chamber(n)*A_star*C_f(n)*Thrust_eff;
    u_e(n) = Thrust(n)/mdot(n);
    if n > 1
        impulse(n) = impulse(n-1) + Thrust(n)*t_step;
    else
        impulse(n) = Thrust(n)*t_step;
    end
    %Update Values in Tank
    m_tank(n+1) = m_tank(n) - mdot(n)*t_step;
    P_tank(n+1) = P_tank(1)*(m_tank(n+1)/m_tank(1))^k;
    T_tank(n+1) = T_tank(1)*(P_tank(n+1)/P_tank(1))^((k-1)/k);
    
    % Update Chamber Temperature
    if P_tank(n+1) > P_c
        P_chamber(n+1) = P_c;
        tstart = 0;
    else
        if tstart == 0
            tstart = time;
            nstart = n;
        end
        P_chamber(n+1) = P_tank(n+1);
    end
    % Update Iteration Values
    n = n+1;
    t(n) = time + t_step;
    time = time + t_step;
end
n = n-1;
figure(1);
plot(t(1:n),P_tank(1:n)/1e6)
title('Tank Total Pressure Variation')
xlabel('Time (s)');
ylabel('Pressure (MPa)');
figure(2);
plot(t(1:n),T_tank(1:n))
title('Tank Total Temperature Variation')
xlabel('Time (s)');
ylabel('Temp (K)');
figure(3);
plot(t(1:n),u_e(1:n))
title('Exit Velocity of Nozzle')
xlabel('Time (s)');
ylabel('Velocity (m/s)');
figure(4)
plot(t(1:n),mdot(1:n))
title('Mass Flow Rate Variation')
xlabel('Time (s)');
ylabel('Mass Flow Rate (kg/s)');
figure(5)
plot(t(1:n),Thrust(1:n))
title('Thrust Variation')
xlabel('Time (s)');
ylabel('Thrust (N)');
figure(6)
plot(t(1:n),P_chamber(1:n)/1e6)
title('Chamber Pressure Variation')
xlabel('Time (s)');
ylabel('Chamber Pressure (MPa)');
figure(7)
plot(t(1:n),C_f(1:n))
title('Thrust Coef Variation')
xlabel('Time (s)');
ylabel('Thrust Coef');
%% Simulate the CGT using only Total Pres Change
% m_expelled_b(1) = 0.0;
% P_tank_b(1) = Tank_total_Pressure_init;
% T_tank_b(1) = Tank_total_Temp_init;
% P_chamber_b(1) = P_c;
% t_b(1) = 0;
% m_tank_b(1) = m_tank_init;
% n = 1;
% t_step = 1;
% time_b = 0;
% 
% while m_expelled_b(n) < 150
%     T_chamber_b(n) = T_tank_b(n);%*(P_chamber(n)/P_tank(n))^((k-1)/k);
%     mdot_b(n) = P_chamber_b(n)*A_star*sqrt(k/R_gas/T_chamber_b(n))*(2/(k+1))^((k+1)/2/(k-1));
%     m_expelled_b(n+1) = m_expelled_b(n) + t_step*mdot_b(n);
%     C_f_b(n) = sqrt(2*k*k/(k-1)*(2/(k+1))^((k+1)/(k-1))*(1 - (P_atm/P_chamber_b(n))^((k-1)/k)));
%     Thrust_b(n) = P_chamber_b(n)*A_star*C_f_b(n)*Thrust_eff;
%     u_e_b(n) = Thrust_b(n)/mdot_b(n);
%     
%     %Update Values in Tank
%     m_tank_b(n+1) = m_tank_b(n) - mdot_b(n)*t_step;
%     T_tank_b(n+1) = T_tank_b(1);
%     P_tank_b(n+1) = m_tank_b(n+1)*R_gas*T_tank_b(n+1)/V_tank;
%     
%     if P_tank_b(n+1) > P_c
%         P_chamber_b(n+1) = P_c;
%     else
%         P_chamber_b(n+1) = P_tank_b(n+1);
%     end
%     % Update Iteration Values
%     n = n+1;
%     t_b(n+1) = time_b + t_step;
%     time_b = time_b + t_step;
% end
%     
% 
% 
% figure(6);
% plot(t_b(1:time_b),P_tank_b(1:time_b)/1e6)
% title('Tank Total Pressure Variation - B')
% xlabel('Time (s)');
% ylabel('Pressure (MPa)');
% figure(7);
% plot(t_b(1:time_b),T_tank_b(1:time_b))
% title('Tank Total Temperature Variation - B')
% xlabel('Time (s)');
% ylabel('Pressure (MPa)');
% figure(8);
% plot(t_b(1:time_b),u_e_b(1:time_b))
% title('Exit Velocity of Nozzle - B')
% xlabel('Time (s)');
% ylabel('Velocity (m/s)');
% figure(9)
% plot(t_b(1:time_b),mdot_b(1:time_b))
% title('Mass Flow Rate Variation - B')
% xlabel('Time (s)');
% ylabel('Mass Flow Rate (kg/s)');
% figure(10)
% plot(t_b(1:time_b),Thrust_b(1:time_b))
% title('Thrust Variation - B')
% xlabel('Time (s)');
% ylabel('Thrust (N)');
% 
% 
% figure(11);
% plot(t(1:time),P_tank(1:time)/1e6,t_b(1:time_b),P_tank_b(1:time_b)/1e6)
% title('Tank Total Pressure Variation - ALL')
% xlabel('Time (s)');
% ylabel('Pressure (MPa)');
% legend('Temp Varies', 'Temp Constant')
% 
% figure(12);
% plot(t(1:time),T_tank(1:time),t_b(1:time_b),T_tank_b(1:time_b))
% title('Tank Total Temperature Variation - ALL')
% xlabel('Time (s)');
% ylabel('Pressure (MPa)');
% legend('Temp Varies', 'Temp Constant')
% 
% figure(13);
% plot(t(1:time),u_e(1:time),t_b(1:time_b),u_e_b(1:time_b))
% title('Exit Velocity of Nozzle - ALL')
% xlabel('Time (s)');
% ylabel('Velocity (m/s)');
% legend('Temp Varies', 'Temp Constant')
% 
% figure(14)
% plot(t(1:time),mdot(1:time),t_b(1:time_b),mdot_b(1:time_b))
% title('Mass Flow Rate Variation - ALL')
% xlabel('Time (s)');
% ylabel('Mass Flow Rate (kg/s)');
% legend('Temp Varies', 'Temp Constant')
% 
% figure(15)
% plot(t(1:time),Thrust(1:time),t_b(1:time_b),Thrust_b(1:time_b))
% title('Thrust Variation - ALL')
% xlabel('Time (s)');
% ylabel('Thrust (N)');
% legend('Temp Varies', 'Temp Constant')

%% OUTPUT
fprintf("Total Time for Propellant to be Expelled\n");
fprintf("Varying Model: %f sec\n",time);
fprintf("Total Mass Expelled: %f kg\n",m_expelled(n-1));
fprintf("Total Impulse: %f Ns\n",impulse(n));
fprintf("Time when BlowDown Starts: %f s\n",tstart);
fprintf("Maximum Thrust: %f N\n",Thrust(nstart));
fprintf("Average Thrust: %f N\n",impulse(n)/time);
fprintf("Liftoff Weight: %f kg  |  %f lbs\n",Thrust(nstart)/5/9.81,Thrust(nstart)/5/9.81*2.2);
