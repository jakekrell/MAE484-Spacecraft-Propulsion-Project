%% Development Information
% MAE 484 Spacecraft Propulsion 
% AAt_to_M_Isentropic_G10.m
% 
% Convert total/local pressure ratio(s) to local Mach number(s).
%
% input[A_At]: local/throat area ratio(s). vector
% input[g]: gamma, ratio of specific heats. scalar
% input[sonic]: "supersonic" or "subsonic" applying to all A_At. string
%
% output[M]: local Mach number(s). vector
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

function [M] = AAt_to_M_Isentropic_G10(A_At,g,sonic)

L = length(A_At);

if strcmpi(sonic,"supersonic")
    M0 = repmat(10,L,1); % initial guess
elseif strcmpi(sonic,"subsonic")
    M0 = repmat(0.5,L,1); % initial guess
else
    disp('Error: input[sonic] must be "supersonic" or "subsonic".')
    return
end
M = zeros(L,1);

gm1 = g-1;
gp1 = g+1;

tol = 1e-10; % tolerance of Newton-Raphson solver
for ii = 1:L
    error = tol+1;

    while abs(error) > tol
        F = (2/gp1*(1+gm1/2*M0(ii)^2))^(gp1/gm1/2) - M0(ii)*A_At(ii); % F = f(M)
        dF = M0(ii)*(2/gp1*(1+gm1/2*M0(ii)^2))^(gp1/2/gm1-1)-A_At(ii); % dF/dM
    
        M(ii) = M0(ii) - F/dF;
        
        error = M(ii) - M0(ii);
        M0(ii) = M(ii);
    end

end

if isrow(A_At)
    M = M'; % transpose output to row based on format of input[A_At]
end

end

