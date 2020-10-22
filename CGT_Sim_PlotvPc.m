%% Setup
clc;
clear;

%% Initial Conditions
thrust = 10; %Newtons
Tank_total_Temp_init = 298.15; %K
Thrust_eff = 0.90; % Efficiency of Nozzle
P_atm = 1.01e5; %Pa - Exit Pressure Equals Atm pressure
% P_c = 2.0684e7; %Pa - Chamber Pressure
k = 1.4; % Specific Heat Ratio
R_gas = 297; %J/kg/K\
C_p = 1.04e3; %J/(kg K) - Heat capacity at a constant pressure

%% Nozzle Exit 
max_iter = 100
P_step = 1e5;
n = 1;
P_c(n) = 2e5;

while n <= max_iter
p_chamber_to_P_exit(n) = P_c(n) / P_atm;

p_exit_to_p_chamber(n) = p_chamber_to_P_exit(n)^-1;

mach_exit(n) = sqrt((2/(k-1)) * ((p_chamber_to_P_exit(n)^((k-1)/k)) - 1));

area_exit_to_throat(n) = (1/mach_exit(n)) * ((1+((k-1)/2)*mach_exit(n)^2)/((k+1)/2))^((k+1)/(2*(k-1)));

temp_exit(n) = Tank_total_Temp_init / (1 + ((k-1)/2) * mach_exit(n)^2);

velocity_exit(n) = sqrt(k*R_gas*temp_exit(n)) * mach_exit(n);
    
area_throat(n) = thrust/(Thrust_eff * P_c(n) * (k*(2/(k+1))^((k+1)/(2*(k-1)))) * sqrt(((2*C_p)/(k*R_gas))*(1 - p_exit_to_p_chamber(n)^((k-1/k)))));

mass_flowrate(n) = (area_throat(n) * P_c(n) /sqrt(k*R_gas*Tank_total_Temp_init)) * (k*(2/(k+1))^((k+1)/(2*(k-1))));

P_c(n+1) = P_c(n) + P_step;

n = n+1;

end
n = n-1;
figure(1)
plot(P_c(1:n),area_exit_to_throat(1:n))
ylabel('Area Ratio')
xlabel('Chamber Pressure (Pa)')
figure(2)
plot(P_c(1:n),mach_exit(1:n))
ylabel('Exit Mach')
xlabel('Chamber Pressure (Pa)')
figure(3)
plot(P_c(1:n),velocity_exit(1:n))
ylabel('Exit Velocity')
xlabel('Chamber Pressure (Pa)')
figure(4)
plot(P_c(1:n),mass_flowrate(1:n))
ylabel('Mass Flowrate (kg/s)')
xlabel('Chamber Pressure (Pa)')

% fprintf('Ae/At = %f \n',area_exit_to_throat)
% fprintf('Throat area = %f m^2\n',area_throat)
% fprintf('Exit area = %f m^2\n',area_throat*area_exit_to_throat)
% fprintf('Exit Mach number = %f \n',mach_exit)
% fprintf('Mass Flowrate = %f kg/s \n',mass_flowrate)
% fprintf('Exit velocity = %f m/s \n',velocity_exit)