p_in = 5.1  ; % [bar]
m_in = 0.039; % [kg/s]
p_out = 1.0 ; % [bar]

Temp = 25; % [C°]
Area = pi*0.0055^2; % [m^2]
gamma = 1.4;

% Density [kg/m^3]
rho_in = (p_in*1e5) / (287.058 * (Temp + 273.15))
%      bar*kg/(m*s^2) / J/(kg*K) * K

% Velocity in x-direction [m/s]
vx_in = m_in/(rho_in*A)
vy_in = 0
%    m^3/(h*bar)

% Energy
E = (p_in*1e5) / ((gamma-1)*rho_in) + 0.5*sqrt(vx_in^2+vy_in^2)