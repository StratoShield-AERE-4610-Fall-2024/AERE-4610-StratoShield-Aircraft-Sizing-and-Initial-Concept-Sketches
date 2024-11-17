clear; clc; clf;
% ========================================================================
%    _______________  ___  __________  _____ __  ______________    ____ 
%   / ___/_  __/ __ \/   |/_  __/ __ \/ ___// / / /  _/ ____/ /   / __ \
%   \__ \ / / / /_/ / /| | / / / / / /\__ \/ /_/ // // __/ / /   / / / /
%  ___/ // / / _, _/ ___ |/ / / /_/ /___/ / __  // // /___/ /___/ /_/ / 
% /____//_/ /_/ |_/_/  |_/_/  \____//____/_/ /_/___/_____/_____/_____/  

% ======================================================================== 
% Initial Sizing | Flight Performance Subteam

% Patricia Belinga (FP Engineer)
% Thomas Housley (FP Lead) 
% Bradley Norwall (CAD Lead)

% ========================================================================

% title made using:
% https://patorjk.com/software/taag/#p=display&f=Slant&t=STRATOSHIELD



%% Reynolds Number Calculations

z = 1500; % [ft] (altitude 150m above Ames from SL reference)

% Standard Atmo values/ratios from Initial Airfoil Sizing source @ 1500 ft
T = 12.2+273.15; % [K]
p_p0 = .947;
rho_rho0 = 0.9568;
mu_mu0 = 0.992;

% sea level conditions for ratios
p0 = 14.7; % [psi]
rho0 = 0.002378; % [slug / ft^3]
mu0 = 3.737e-7; % [slug / ft * s]
Cruise.Vel = 58.6667; % [ft / s] = 40 mph 
Cruise.Vel_SI = 17.8816; % [m / s]

% chord was used for Reynold's, span is just under max wingspan and kinda
% = 10 max Aspect Ratio recommendation from class THESE CHANGE BASED ON
% SIZING PARAMETERS FROM STRUCTURES TEAM *******************************
% IF THAT HAPPENS THEN NEW XFLR5 DATA NEEDS TO BE AQUIRED
chord = 8 / 12; % [ft]
span = 85 / 12; % [ft]

% now getting actual values 
rho = rho_rho0 * rho0;
mu = mu_mu0 * mu0;

% Reynolds # at cruise speed (40mph)
Cruise.Re = (rho * Cruise.Vel * chord) / mu;

% using SI because mach is unitless, verified that it matches cruise in MPH
a = sqrt(T * 1.4 * 287);
Cruise.Mach = Cruise.Vel_SI / a;

%% Loading CSV XFLR5 data

% ========================================================================
% IF YOU ARE RUNNING THIS ON YOUR OWN COMPUTER, THE XFLR5 CSV FILES WILL
% PROBABLY HAVE A DIFFERENT PATH. IF THEY HAVE THE SAME FILENAME AS
% FORMATTED BELOW:

% -RIGHT CLICK ON THE FILE IN YOUR FILE EXPLORER
% -SELECT "COPY AS PATH" (CTRL+SHIFT+C)
% -PASTE PATH INTO directory VARIABLE BELOW 
% -REMOVE THE FILENAME FROM THE PATH STRING YOU PASTED 

% then you should be able to use this code
% ========================================================================

naca4412_filename = "NACA_4412_data.csv";
naca6409_filename = "NACA_6409_data.csv";
mh60_filename = "MH60_data.csv";
e423_filename = "E423_data.csv";

directory = "C:\Skewl\AerE4610\Airfoil XFLR5 initial sizing data\";

% creating path string for each airfoil

naca4412_path = strjoin(directory + naca4412_filename);
naca6409_path = strjoin(directory + naca6409_filename);
mh60_path = strjoin(directory + mh60_filename);
e423_path = strjoin(directory + e423_filename);

% loading variable columns into table so its easy to identify columns

naca4412_data = readtable(naca4412_path, ...
    "Delimiter", ",", "VariableNamingRule", "preserve");
naca6409_data = readtable(naca6409_path, ...
    "Delimiter", ",", "VariableNamingRule", "preserve");
mh60_data = readtable(mh60_path, ...
    "Delimiter", ",", "VariableNamingRule", "preserve");
e423_data = readtable(e423_path, ...
    "Delimiter", ",", "VariableNamingRule", "preserve");

% this is for loops later on dw
n_naca4412 = size(naca4412_data);
n_naca6409 = size(naca6409_data);
n_mh60 = size(mh60_data);
n_e423 = size(e423_data);


% recreating XFLR5 plots for viewing
% this can be done for tail airfoils too so I can add that code

figure(1)
plot(naca4412_data.CD, naca4412_data.CL)
hold on
plot(naca6409_data.CD, naca6409_data.CL)
plot(mh60_data.CD, mh60_data.CL)
% plot(e423_data.CD, e423_data.CL)
grid on
title("C_l vs. C_d")
legend("NACA 4412", "NACA 6409", "MH60")%,"E423 h-stabilizer")
xlabel("C_d")
ylabel("C_l")
hold off

% time for curve fitting (throwing up in my mouth rn)
% im sorry its hardcoded on the rows and columns
fit_naca4412 = polyfit( ...
    naca4412_data.alpha(5:22), naca4412_data.CL(5:22), 1);
fit_naca6409 = polyfit( ...
    naca6409_data.alpha(5:22), naca6409_data.CL(5:22), 1);
fit_mh60 = polyfit(mh60_data.alpha(5:22), mh60_data.CL(5:22), 1);
fit_e423 = polyfit(e423_data.alpha(5:22), e423_data.CL(5:22), 1);

figure(2)
plot(naca4412_data.alpha, naca4412_data.CL)
hold on
plot(naca6409_data.alpha, naca6409_data.CL)
plot(mh60_data.alpha, mh60_data.CL)
% plot(e423_data.alpha, e423_data.CL)

% plotting fits (throwing up in my mouth again)
plot(naca4412_data.alpha, ...
    fit_naca4412(1) .* naca4412_data.alpha + fit_naca4412(2), ...
    "LineStyle", "--", "Color", "black")
plot(naca6409_data.alpha, ...
    fit_naca6409(1) .* naca6409_data.alpha + fit_naca6409(2), ...
    "LineStyle", "--", "Color", "black")
plot(mh60_data.alpha, ...
    fit_mh60(1) .* mh60_data.alpha + fit_mh60(2), ...
    "LineStyle", "--", "Color", "black")
% plot(e423_data.alpha, ...
%     fit_e423(1) .* e423_data.alpha + fit_e423(2), ...
%     "LineStyle", "--", "Color", "black")

grid on
title("AoA vs. C_l")
xlabel("AoA (degrees)")
ylabel("C_l")
xline(0)
yline(0)
legend("NACA 4412", "NACA 6409", "MH60") % , "E423 h-stabilizer")
hold off

figure(3)
plot(naca4412_data.alpha, naca4412_data.Cm)
hold on
plot(naca6409_data.alpha, naca6409_data.Cm)
plot(mh60_data.alpha, mh60_data.Cm)
% plot(e423_data.alpha, e423_data.Cm)
grid on
title("AoA vs. C_m")
xlabel("AoA (degrees)")
ylabel("C_m")
xline(0)
yline(0)
xlim([-4,20])
legend("NACA 4412", "NACA 6409", "MH60") % , "E423 h-stabilizer")
hold off

% ok ze curves hath been createth
% i only did curves again so i could have the slopes for later hehe
% from the polyfits can then identify L/D 

%% LIFT/DRAG CALCULATIONS (OPTIMIZATION)
% all of the (important) calculated data thats importnant is dumped into
% calc struct

% Cl/Alpha 2D slopes
calc.naca4412.Lalpha = fit_naca4412(1);
calc.naca6409.Lalpha = fit_naca6409(1);
calc.mh60.Lalpha = fit_mh60(1);
calc.e423.Lalpha = fit_e423(1);

% Cl^(3/2)/C_d ratio matrices
% aere 261 notes on max battery prop endurance (see teams folder w slides)
% calc.naca4412.LD_aloft_values = (naca4412_data.CL).^(3 / 2) ...
%     ./ naca4412_data.CD;
% calc.naca6409.LD_aloft_values = (naca6409_data.CL).^(3 / 2) ...
%     ./ naca6409_data.CD;
% calc.mh60.LD_aloft_values = (mh60_data.CL).^(3 / 2) ./ mh60_data.CD;
% calc.e423.LD_aloft_values = (e423_data.CL).^(3 / 2) ./ e423_data.CD;

% maximum Cl^(3/2)/C_d ratios
% calc.naca4412.LDmax_aloft = max(calc.naca4412.LD_aloft_values);
% calc.naca6409.LDmax_aloft = max(calc.naca6409.LD_aloft_values);
% calc.mh60.LDmax_aloft = max(calc.mh60.LD_aloft_values);
% calc.e423.LDmax_aloft = max(calc.e423.LD_aloft_values);

% C_l maxes and the row location of the max to later find corresponding
% stall angle 
[calc.naca4412.CLmax, calc.naca4412.stall_index]= max(naca4412_data.CL);
[calc.naca6409.CLmax, calc.naca6409.stall_index] = max(naca6409_data.CL);
[calc.mh60.CLmax, calc.mh60.stall_index]= max(mh60_data.CL);
[calc.e423.CLmax, calc.e423.stall_index]= max(e423_data.CL);

% C_d maxes now 
calc.naca4412.CDmax = max(naca4412_data.CD);
calc.naca6409.CDmax = max(naca6409_data.CD);
calc.mh60.CDmax = max(mh60_data.CD);
calc.e423.CDmax = max(e423_data.CD);


calc.naca4412.LD_aloft_max = (calc.naca4412.CLmax)^(3 / 2) ...
    / calc.naca4412.CDmax;
calc.naca6409.LD_aloft_max = (calc.naca6409.CLmax)^(3 / 2) ...
    / calc.naca6409.CDmax;
calc.mh60.LD_aloft_max = (calc.mh60.CLmax)^(3 / 2) / calc.mh60.CDmax;
calc.e423.LD_aloft_max = (calc.e423.CLmax)^(3 / 2) / calc.e423.CDmax;

% afformentioned alpha_stall
calc.naca4412.alpha_stall = naca4412_data.alpha(calc.naca4412.stall_index);
calc.naca6409.alpha_stall = naca6409_data.alpha(calc.naca6409.stall_index);
calc.mh60.alpha_stall = mh60_data.alpha(calc.mh60.stall_index);
calc.e423.alpha_stall = e423_data.alpha(calc.e423.stall_index);

% HOWEVER, THESE RESULTS ARE BASED ON CL_MAX
% can also get "linear stall alpha"
% this value is lower than other derived alpha_stall but may be good
% comparison?
% best way i can think of deriving it in a repeatable way is having a
% tolerance between linear fit line and actual XFLR5 data

tolerance = 0.1;

difference_naca4412 = ...
    (fit_naca4412(1) .* naca4412_data.alpha + fit_naca4412(2)) ...
    - naca4412_data.CL;
difference_naca6409 = ...
    (fit_naca6409(1) .* naca6409_data.alpha + fit_naca6409(2)) ...
    - naca6409_data.CL;
difference_mh60 = (fit_mh60(1) .* mh60_data.alpha + fit_mh60(2)) ...
    - mh60_data.CL;
difference_e423 = (fit_e423(1) .* e423_data.alpha + fit_e423(2)) ...
    - e423_data.CL;

calc.naca4412.linearstall = 0;
calc.naca6409.linearstall = 0;
calc.mh60.linearstall = 0;
calc.e423.linearstall = 0;
% gross but short nested loops sorry this code is O(n^100)
for i = 1:n_naca4412(1)
    if difference_naca4412(i) >  tolerance
        calc.naca4412.linearstall = naca4412_data.alpha(i);
        break;
    end
end
for i = 1:n_naca6409(1)
    if difference_naca6409(i) > tolerance
        calc.naca6409.linearstall = naca6409_data.alpha(i);
        break;
    end
end
for i = 1:n_mh60(1)
    if difference_mh60(i) > tolerance
        calc.mh60.linearstall = mh60_data.alpha(i);
        break;
    end
end
for i = 1:n_e423(1)
    if difference_e423(i) > tolerance
        calc.e423.linearstall = e423_data.alpha(i);
        break;
    end
end

% ok so basically after doing all that they are almost the same but having
% 2 is good ig? idk i am tired

% GLORIOUS DATA AQUIRED
% time for consideration

% i would start controls but need moments of inertia and I will throw up
% making 4 billion assumptions for that so no


disp("NACA 4412:"); disp(calc.naca4412)
disp("NACA 6409:"); disp(calc.naca6409)
disp("MH60:"); disp(calc.mh60)
disp("E423:"); disp(calc.e423)


%% Reynolds Number Range (RUN CLEAR, CLC BEFORE RUNNING)

% Need a range of Re values for flight at 30mph - 60mph on ground and
% cruise altitude

ground.z = 304; % [m] (~altitude on ground @ ames from SL) (1000ft)
cruise.z = 457; % [m] (altitude ~150m above ames from SL) (1500ft)

% Standard Atmo values & ratios from Initial Airfoil Sizing source @ cruise
cruise.T = 12.2 + 273.15; % [K]
cruise.p_p0 = 0.947;
cruise.rho_rho0 = 0.9568;
cruise.mu_mu0 = 0.992;

% Standard Atmo values & ratios from Initial Airfoil Sizing source @ ground
ground.T = 13.2 + 273.15; % [K]
ground.p_p0 = 0.9644;
ground.rho_rho0 = 0.9711;
ground.mu_mu0 = 0.9947;

% sea level conditions for ratios
p0 = 101325; % [Pa]
rho0 = 1.225; % [kg / m^3]
mu0 = 1.789e-5; % [kg / (m * s)]


vel = [30, 40, 50, 60] ./ 2.237; % mph --> m/s conversion 

% chord was used for Reynold's, span is just under max wingspan and about
% equal 10 max Aspect Ratio recommendation from Class 

% horizontal stabilizer chord was arbitrary choice and may be subject to
% change

chord_wing = 8 / 39.37; % in --> meters
span_wing = 85 / 39.37; % in --> meters

chord_hStab = 6 / 39.37; % in --> meters
span_hor = 15 / 39.37; % in --> meters


% now getting actual values
% cruise conditions
cruise.rho = cruise.rho_rho0 * rho0;
cruise.mu = cruise.mu_mu0 * mu0;

% ground conditions
ground.rho = ground.rho_rho0 * rho0;
ground.mu = ground.mu_mu0 * mu0;

% Reynolds #s
cruise.Re_wing = (cruise.rho .* vel .* chord_wing) ./ cruise.mu;
cruise.Re_stab = (cruise.rho .* vel .* chord_hStab) ./ cruise.mu;
ground.Re_wing = (ground.rho .* vel .* chord_wing) ./ ground.mu;

% finding machs
cruise.a = sqrt(cruise.T * 1.4 * 287);
cruise.mach = vel ./ cruise.a;
ground.a = sqrt(ground.T * 1.4 * 287);
ground.mach = vel ./ ground.a;

disp("Ground Conditions to be input to XFLR5")
disp(ground)
disp("Cruise Conditions to be input to XFLR5")
disp(cruise)
