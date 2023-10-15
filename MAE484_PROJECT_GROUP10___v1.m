%% Development Information
% MAE 484 Spacecraft Propulsion - Midterm Project
% 
% Calculate various items throughout the project tasks.
% 
% 
% Developers Contact Information:
% Vinicius D. Ferreira, Jacob P. Krell
% Aerospace Engineering Undergraduate Students
% Statler College of Engineering & Mineral Resources
% Dept. Mechanical and Aerospace Engineering
% West Virginia University (WVU)
% vdf00001@mix.wvu.edu, jpk0024@mix.wvu.edu
%
%
%
% Development History
% Date              Developer(s)     Comments
% ---------------   -------------    --------------------------------
% Sept. 29, 2023    VDF, JPK         Setting up givens and tasks 2-4.
% Sept. 30, 2023    VDF, JPK         Completion of tasks 2-3.
%
%% Initialization of Given Data

func = '\func';
folder = pwd;
if folder(end-4:end) ~= func
    addpath(append(pwd,func)) % add folder with functions to path
end

clearvars, clc, close all

% Nozzle and Propellant (NTO/A50) Characteristics:
Pc = 8.15 * 101325; % Pa, chamber pressure (actual)
Tc_ideal = 3151; % K, ideal chamber temperature
F_ideal = 16*1e3; % N, ideal thrust
Isp_ideal = 339; % s, ideal specific impulse
m_inert = 2150; % kg, inert mass
O_F = 1.6/1;
m_pay = 245; % kg, payload mass
dv = 2*1e3; % m/s, change in velocity (delta-v)
E = 46; % area ratio

% Correction Factors:
n_cstar = 0.95;
n_cF = 0.96;
n_v = n_cF * n_cstar;

% Gas Characteristics:
R = 408.6; % J/(kg*K)
g = 1.14; % ratio of specific heats
gp1 = g + 1;
gm1 = g - 1;

% Gravity:
g0 = 9.8066; % m/s^2

% -------------------------------Task 2------------------------------------
% Gross Vehicle Sizing

Isp = n_v * Isp_ideal; % s, actual specific impulse
c_ideal = Isp * g0; % m/s, effective exhaust velocity
c = Isp * g0; % m/s, effective exhaust velocity

m_f = m_inert + m_pay; % kg, final mass
m_prop = m_f * (exp(dv/c_ideal) - 1); % kg, propellant mass
m_i = m_f + m_prop; % kg, initial mass
m_t = m_pay + m_prop + m_inert; % kg, total mass
r_pay = m_pay / m_t; % unitless, payload mass ratio
zeta_inert = m_inert / (m_inert + m_prop); % unitless, Inert mass fraction

I_t = Isp * m_prop * g0; % N * s, Total Impulse
m_dot_ideal = F_ideal / c_ideal; % kg / s, mass flow rate ideal
t_burn = I_t / F_ideal; % s, burn duration


% -------------------------------Task 3------------------------------------
% Isentropic Nozzle Flow

% Part (a):

Tc_Tt = M_to_T0T_Isentropic_G10(1,g);
Tt_ideal = Tc_ideal/Tc_Tt; % K, ideal throat temperature

Pc_Pt = M_to_P0P_Isentropic_G10(1,g); % ... see below

Tt = n_cstar^2*Tt_ideal; % K, throat temperature (actual)
% Pt = ; % Pa, throat pressure (actual) ... see below
vt = sqrt(g*R*Tt); % m/s, throat velocity (actual)

% Part (b):

Me = AAt_to_M_Isentropic_G10(E,g,'supersonic'); % exit Mach number

Tc_Te = M_to_T0T_Isentropic_G10(Me,g);
Te_ideal = Tc_ideal/Tc_Te; % K, ideal exit temperature

Pc_Pe = M_to_P0P_Isentropic_G10(Me,g); % ... see below

Te = n_cstar^2*Te_ideal; % K, exit temperature (actual)
% Pe = ; % Pa, exit pressure (actual) ... see below
ve = sqrt(g*R*Te); % m/s, exit velocity (actual)

% Part (c):

cstar_ideal = cstar_Isentropic_G10(g,R,Tc_ideal); % m/s, ideal characteristic velocity

cstar = n_cstar*cstar_ideal; % m/s, characteristic velocity (actual)
cF = c/cstar; % thrust coefficient (actual)

% Part (d):

mdot = F_ideal / c; % kg/s, mass flow rate (actual)

At = cstar*mdot/Pc; % m^2, throat area
Ae = E*At; % m^2, exit area

dt = 2*sqrt(At/pi); % m, throat diameter
de = 2*sqrt(Ae/pi); % m, exit diameter

% Returning to P_actual calculations with known At:

mdot_ideal = F_ideal/c_ideal;
Pc_ideal = cstar_ideal*mdot_ideal/At; % NOTE Pc_ideal=Pc MIGHT BE TRUE BY DEFINITION !!!!!!!!!!!!!!!!

Pt_ideal = Pc_ideal/Pc_Pt; % Pa, ideal throat pressure
% Pt = ; % Pa, throat pressure (actual)

Pe_ideal = Pc_ideal/Pc_Pe; % Pa, ideal exit pressure
% Pe = ; % Pa, exit pressure (actual)


% -------------------------------Task 4------------------------------------
% Propellant Tank Sizing

% Two tanks, both spherical
% Oxidizer - NTO (Nitrogen Tetroxide), Propellant - A50 (Aerozine 50/50)
rho_ox = 1440;     % kg/m^3, NTO density
rho_fuel = 884.55; % kg/m^3, A50 density

n_vol = 0.95; % tank expulsion efficiency (for each tank)

% Pressure loss
v_ox = 16.5;      % m/s, oxidizer fluid flow 
v_fuel = 21.07;   % m/s, fuel fluid flow
P_feed = 40000;   % Pa, line pressure loss
P_inj = 0.3 * Pc; % Pa, injector pressure loss

P_dyn_ox = (rho_ox * (v_ox^2)) / 2; % Pa, Dynamic pressure loss oxidizer
P_dyn_fuel = (rho_fuel * (v_fuel^2)) / 2; % Pa, Dynamic pressure loss fuel

P_He = 25e6; % Pa, Pressure fed (Helium)
rho_He = 0.1786; % kg/m^3, Helium density
g_He = 1.66; % Ratio of specific heat Helium

fuel_m = m_prop / (O_F + 1); % kg, fuel mass
ox_m = (m_prop * O_F) / (O_F + 1); % kg, oxidizer mass

fuelTank_V = fuel_m / rho_fuel; % m^3, fuel tank volume
oxTank_V  = ox_m / rho_ox; % m^3, oxidizer tank volume

oxTank_r = ((3 * oxTank_V) / (4 * pi)) ^ (1/3); % m, radius of oxidizer tank
fuelTank_r = ((3 * fuelTank_V) / (4 * pi)) ^ (1/3); % m, radius of fuel tank

oxTank_P = Pc + P_dyn_ox + P_feed + P_inj; % Pa, Oxidizer tank pressure
fuelTank_P = Pc + P_dyn_fuel + P_feed + P_inj; % Pa, Fuel tank pressure

Ti_Tf_ratio = P_to_T_Isentropic_G10(P_He/fuelTank_P,g);

% m^3, Pressurant gas tank volume
totalTank_V = (fuelTank_P*fuelTank_V)/(P_He*Ti_Tf_ratio - fuelTank_P); 
% kg, Pressurant gas tank mas
totalTank_m = totalTank_V*rho_He; 

% -------------------------------Task 5------------------------------------
% Combustion Chamber Sizing

Lstar = 0.7; % m, characteristic length
Mc = 0.1654; % chamber Mach
 
% Part (a)

% contraction ratio
Ec = (1/Mc)*sqrt( ((2/gp1)*(1 + (gm1/2)*Mc^2))^(gp1/gm1) ); 
 
% Part (b)

dc = sqrt(4*Ec*At/pi); % m, Chamber Diameter
rc = dc/2; % m, Chamber Radius
Ac = (pi*dc^2)/4; % m^2, Chamber Area
 
% Part ©

Lc = Lstar/Ec; % m, Chamber Length
Vc = At*Lstar; % m^3, Chamber Volume

Lc_dc = Lc/dc; % Good, fall into range 0.5 < Lc/dc < 2.5

% -------------------------------Task 6------------------------------------
% Conical Nozzle Design

alp = 15; % deg, half angle for a conical nozzle

% Part a
Ln = (de - dt)/(2*tand(alp)); % m, Nozzle Length
lambda = ((1+cosd(alp))/2); % percent, angle correction factor

% Part b
rt = dt/2; % m, throat radius 
rtd = 0.382*rt; % m, throat outlet radius
Ln_r = (rt*(sqrt(E)-1) + rtd*(secd(alp)-1)) / (tand(alp));

% Part c

% plt == plot controls = . . .
% . . . [plt_global,plt_sect,plt_mid,plt_radii,plt_nTyp,plt_nPnt,plt_sTyp,
% . . . plt_sPnt,plt_mTyp,plt_mPnt,plt_rTyp,plt_rPnt,plt_buffer]
% plt = {1,1,1,1,'k-',2,'b--',1,'k-.',0.5,'r-.',0.5,0.05}; % v1
plt = {1,1,1,0,'k-',4,'k-.',0.5,'k-.',0.5,'r-.',0.5,0.05}; % v2
% Nozzle = CreateNozzleConical2D_G10(rc,Lc,1.5*rt,rt,rtd,alp,Ln_r,1e2,plt,0);
Nozzle = CreateNozzleConical2D_StraightConvergent_G10(rc,Lc,1.5*rt,135,rt,rtd,alp,Ln_r,1e2,plt,0);

% Part d
thi = 47;
thf = 10;
