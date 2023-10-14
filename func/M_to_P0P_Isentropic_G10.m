%% Development Information
% MAE 484 Spacecraft Propulsion 
% M_to_P0P_Isentropic_G10.m
% 
% Convert local Mach number(s) to total/local pressure ratio(s).
%
% input[M]: local Mach number(s). vector
% input[g]: gamma, ratio of specific heats. scalar
%
% output[P0_P]: total/local pressure ratio(s). vector
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
% Sept. 12, 2023    J. Krell         For nozzle isentropic eq.'s catalog
%
%%

function [P0_P] = M_to_P0P_Isentropic_G10(M,g)

gm1 = g-1;

P0_P = (1+gm1/2*M.^2).^(g/gm1);


end

