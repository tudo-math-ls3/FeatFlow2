% Pressure and mass flow rate at the inlet
p_in = 510000 ; % [Pa]
m_in = 0.039; % [kg/s]

% Area at the inlet
Area = pi*0.0055^2; % [m^2]

% Pressure at the outlet
p_0  = 100000 ; % [Pa]

% Temperature and gas constant
Temp = 25; % [C°]
gamma = 1.4;

% --- Physical quantities ---
disp('Physical quantities')

% Density [kg/m^3]
rho_in = (p_in) / (287.058 * (Temp + 273.15))
rho_0  = (p_0) / (287.058 * (Temp + 273.15))

% Velocity in x-direction [m/s]
vx_in = m_in/(rho_in*Area)
vy_in = 0

% Pressure
p_in
p_0

% Energy
E_in = (p_in) / (gamma-1) + 0.5*rho_in*(vx_in^2+vy_in^2)
E_0  = (p_0) / (gamma-1)


% Reference values
rho_r = rho_in
vel_r = vx_in
len_r = 1.0


% --- Dimensionless quantities ---
disp('Dimensionless quantities')

% Density
rho_in = rho_in/rho_r
rho_0  = rho_0/rho_r

% Velocity
vx_in = vx_in/vel_r
vy_in = vy_in/vel_r

% Pressure
p_in = p_in/(rho_r*vel_r^2)
p_0  = p_0/(rho_r*vel_r^2)

% Energy
E_in = E_in/(rho_r*vel_r^2)
E_0  = E_0/(rho_r*vel_r^2)