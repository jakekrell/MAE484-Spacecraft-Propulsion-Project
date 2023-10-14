%% Development Information
% MAE 484 Spacecraft Propulsion 
% cstar_Isentropic_G10.m
% 
% Convert total temperature(s) (K) to characteristic velocity(s).
%
% input[g]: gamma, ratio of specific heats. scalar
% input[R]: gas constant of propellant. scalar
% input[T0]: total temperature(s) (K). vector
%
% output[cstar]: characteristic velocity(s). vector
% 
% Assumptions:
% (1) Isentropic expansion
% (2) Flow within 1D nozzle
% 
% 
% Primary Developer Contact Information:
% Jacob P. Krell [Project Group 10]
% Aerospace Engineering Undergraduate Student
% Statler College of Engineering & Mineral Resources
% Dept. Mechanical and Aerospace Engineering
% West Virginia University (WVU)
% jpk0024@mix.wvu.edu
%
%
%
% Development History
% Date              Developer        Comments
% ---------------   -------------    --------------------------------
% Sept. 21, 2023    J. Krell         For nozzle isentropic eq.'s catalog
%
%%

function [cstar] = cstar_Isentropic_G10(g,R,T0)

gm1 = g-1;
gp1 = g+1;

cstar = sqrt(g*R*T0)/g/sqrt((2/gp1)^(gp1/gm1));


end

