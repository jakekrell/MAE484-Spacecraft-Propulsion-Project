%% Development Information
% MAE 484 Spacecraft Propulsion 
% P_to_T_Isentropic_Ferreira.m
% 
% Convert Pressure ratio to a Temperature ratio at any two arbitrary
% locations in the nozzle
%
% input[Px_Py]: Pressure Ratio
% input[gama] : Ratio of specific heats
%
% output[Tx_Ty]: Temperature Ratio based on Pressure ratio given
% 
% Assumptions:
% (1) Isentropic Flow
% 
% 
% Primary Developer Contact Information:
% Vinicius Dunker Ferreira, Aerospace Engineering
% Undergraduate Student
% Statler College of Engineering & Mineral Resources
% Dept. Mechanical and Aerospace Engineering
% West Virginia University (WVU)
% vdf00001@mail.wvu.edu
%
%
% Development History
% Date              Developer            Comments
% ---------------   -------------        --------------------------------
% Sept. 04, 2023    Vinicius Ferreira    Initial implemention
%
%%

function [Tx_Ty] = P_to_T_Isentropic_Ferreira(Px_Py, gama)

    gm1 = gama - 1;
    Tx_Ty = Px_Py ^ ( gm1 / gama );

end