%% Setup
clc;
clear;

%% Initial Conditions
thrust = 10; %Newtons
Tank_total_Temp_init = 298.15; %K
Thrust_eff = 0.90; % Efficiency of Nozzle
P_atm = 1.01e5; %Pa - Exit Pressure Equals Atm pressure
P_c = 4e7; %Pa - Chamber Pressure
k = 1.4; % Specific Heat Ratio
R_gas = 297; %J/kg/K\
C_p = 1.04e3; %J/(kg K) - Heat capacity at a constant pressure

%% Nozzle Exit 

p_chamber_to_P_exit = P_c / P_atm;
p_exit_to_p_chamber = p_chamber_to_P_exit^-1;

% solving for exic mach with pressure ratio
mach_exit = sqrt((2/(k-1)) * ((p_chamber_to_P_exit^((k-1)/k)) - 1));

% Finding exit to throat area ratio using mach exit
area_exit_to_throat = (1/mach_exit) * ((1+((k-1)/2)*mach_exit^2) ...
    /((k+1)/2))^((k+1)/(2*(k-1)));

% Solving for exit temperature to get the exit velocity
temp_exit = Tank_total_Temp_init / (1 + ((k-1)/2) * mach_exit^2);
velocity_exit = sqrt(k*R_gas*temp_exit) * mach_exit;

% Solving for throat area
area_throat = thrust/(Thrust_eff * P_c * (k*(2/(k+1))^((k+1)/(2*(k-1))))...
    * sqrt(((2*C_p)/(k*R_gas))*(1 - p_exit_to_p_chamber^((k-1/k)))));
% area_throat = thrust/(Thrust_eff * P_c * (k*(2/(k+1))^((k+1)/(2*(k-1))))...
%      * velocity_exit);
% 
% 
% Using throat area to find mass flowrate
mass_flowrate = (area_throat * P_c /sqrt(k*R_gas*Tank_total_Temp_init)) ...
    * (k*(2/(k+1))^((k+1)/(2*(k-1))));
% v_exit_2 = sqrt(((2*C_p*Tank_total_Temp_init))*(1 - p_exit_to_p_chamber^((k-1/k))))

% Printing values to command window
fprintf('Ae/At = %f \n',area_exit_to_throat)
fprintf('Throat area = %f m^2\n',area_throat)
fprintf('Exit area = %f m^2\n',area_throat*area_exit_to_throat)
fprintf('Exit Mach number = %f \n',mach_exit)
fprintf('Mass Flowrate = %f kg/s \n',mass_flowrate)
fprintf('Exit velocity = %f m/s \n',velocity_exit)



