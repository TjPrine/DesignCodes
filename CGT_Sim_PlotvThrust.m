%% Setup
clc;
clear;

%% Initial Conditions
%thrust = 200; %Newtons
Tank_total_Temp_init = 298.15; %K
Thrust_eff = 0.90; % Efficiency of Nozzle
P_atm = 1.01e5; %Pa - Exit Pressure Equals Atm pressure
P_c = 2.0684e7; %Pa - Chamber Pressure
k = 1.4; % Specific Heat Ratio
R_gas = 297; %J/kg/K\
C_p = 1.04e3; %J/(kg K) - Heat capacity at a constant pressure

%% Nozzle Exit 
max_iter = 100
n = 1;
thrust_init = 0;
thrust_step = 1;
thrust(n) = thrust_init;

p_chamber_to_P_exit = P_c / P_atm;
p_exit_to_p_chamber = p_chamber_to_P_exit^-1;
mach_exit = sqrt((2/(k-1)) * ((p_chamber_to_P_exit^((k-1)/k)) - 1));

area_exit_to_throat = (1/mach_exit) * ((1+((k-1)/2)*mach_exit^2)/((k+1)/2))^((k+1)/(2*(k-1)));

temp_exit = Tank_total_Temp_init / (1 + ((k-1)/2) * mach_exit^2);

velocity_exit = sqrt(k*R_gas*temp_exit) * mach_exit;

while n <= max_iter
    
area_throat(n) = thrust(n)/(Thrust_eff * P_c * (k*(2/(k+1))^((k+1)/(2*(k-1)))) * sqrt(((2*C_p)/(k*R_gas))*(1 - p_exit_to_p_chamber^((k-1/k)))));

mass_flowrate = (area_throat * P_c /sqrt(k*R_gas*Tank_total_Temp_init)) * (k*(2/(k+1))^((k+1)/(2*(k-1))));

thrust(n+1) = thrust(n) + thrust_step;

n = n+1;

end
n = n-1;
figure(1)
plot(thrust(1:n),area_throat(1:n))
ylabel('Throat Area')
xlabel('Thrust')

fprintf('Ae/At = %f \n',area_exit_to_throat)
fprintf('Throat area = %f m^2\n',area_throat)
fprintf('Exit area = %f m^2\n',area_throat*area_exit_to_throat)
fprintf('Exit Mach number = %f \n',mach_exit)
fprintf('Mass Flowrate = %f kg/s \n',mass_flowrate)
fprintf('Exit velocity = %f m/s \n',velocity_exit)