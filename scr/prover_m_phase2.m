function [varargout] = prover_m_phase2(Cloud_Input,Hopper,Sediment,Para,d,rho_amb,istep,qe,amb,z,dt,varargin)
%% prover_m_phase2: Dynamic Collapse
% Made by:          Jannek Gundlach
% Date of Change:   13.01.2023
% License:          GNU GPL
% Definition:       Sediment spreading radially on the ground after impact.
%                   Calculating the change in x, y, and z direction based on the conservation of 
%                   energy per Runge-Kutta 4 (RK4) fraction of the time step 
% Phase 2 
% End:              The Cloud becomes passive, ambient velocity is larger than the motion of the cloud (checked in main)

%   Define input parameter
defaultiz = 1;
defaultg = 9.81;


p = inputParser;
p.FunctionName = 'prover_m_phase2';
p.CaseSensitive = true;
p.StructExpand = true;
validStruct = @(x) isstruct(x);
validDouble = @(x) isnumeric(x) && isfloat(x);
validPosDouble = @(x) isnumeric(x) && isfloat(x) && all(x>=0);
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
validMatrix = @(x) ismatrix(x) && all(x>=0);
addRequired(p,'Cloud_Input',validDouble);
addRequired(p,'Hopper',validStruct);
addRequired(p,'Sediment',validStruct);
addRequired(p,'Para',validStruct);
addRequired(p,'d',validPosDouble);
addRequired(p,'rho_amb',validPosDouble);

addRequired(p,'qe',validDouble);
addRequired(p,'amb',validDouble);
addRequired(p,'istep',validScalarPosNum);
addRequired(p,'z',validMatrix);
addRequired(p,'dt',validPosDouble);
addOptional(p,'iz',defaultiz,validPosDouble);
addOptional(p,'g',defaultg,validPosDouble);

parse(p,Cloud_Input,Hopper,Sediment,Para,d,rho_amb,qe,amb,istep,z,dt,varargin{:});

k = length(Sediment.name(:));
%% Input preparation
% Check for valid input
if Cloud_Input(3) > d
    Cloud_Input(3) = d;
elseif Cloud_Input(3) == 0.0
    'Input error,please check your settings'
else
    Cloud_Input(3) = Cloud_Input(3);
end
% iterate for the right depth coordinate and defined ambient density
iz = p.Results.iz;
while Cloud_Input(3) > z(iz)
    iz = iz + 1;
end

% Interpolate ambient density at cloud depth
rho_amb_tmp = rho_amb(iz) + (Cloud_Input(3) - z(iz)) * (rho_amb(2) - rho_amb(1)) / (z(2) - z(1));

%% assign variables to easy-use variables

q_e = qe/dt;                        % Entrainment stream (volume per time step)
a = Cloud_Input(7);                 % Cloud(timestep-1).a
b = Cloud_Input(8);                 % Cloud(timestep-1).b
c = Cloud_Input(9);                 % Cloud(timestep-1).c
volumen = Cloud_Input(10);


%% Calculation of the Dynamic Collapse Output per time step

rho_bulk = Cloud_Input(4) / volumen;                                                % Cloud density

% Terms for the roots-equation to determine the spreading in x-direction
roots_term1 = (1 + (b^2 / c^2) + (a^2 / b^2) + (2*a^2 / c^2) + (a^2 * b^2 / c^4));  % Term 1 for roots
roots_term2 = 2 * q_e / volumen * (1/b + b/c^2);                                    % Term 2 for roots
roots_term3 = 10 * Cloud_Input(5) / (rho_bulk * volumen) - a^2 * q_e^2 / volumen^2; % Term 3 for roots
[tmp] = roots([roots_term1 -roots_term2 -roots_term3]);                             % Applying roots with two results
if tmp(2) <=0                                                                       % Check if one of the results is negative or zero
    DynColl_Output(8) = tmp(1);                                                     % If true than the other result is the wanted one
else
    [val,idx] = min(abs(tmp-1));                                                    % If false the result closer to one is the wanted value
    DynColl_Output(8) = tmp(idx);                                                   % Take that value as "db / dt"
end

% Determine the spreading value for y and z
DynColl_Output(9) = (b/c) * DynColl_Output(8);                                      % "dc / dt"
DynColl_Output(7) = a * (q_e/volumen - DynColl_Output(8)/b - DynColl_Output(9)/c);  % "da / dt" the height will reduce

% Calculate the work necessary in the next time step to get further spreading
d_bottom = (DynColl_Output(8)^2 + DynColl_Output(9)^2) * sqrt(DynColl_Output(8)^2 + DynColl_Output(9)^2);
% Work components for friction based on the ground interface
W_f = 0.25 * rho_bulk * Para.CFRIC * pi * b * c *(d_bottom);                        % Friction part of the work
% Work components for drag based on the cloud surface
W_d = 0.25 * rho_amb_tmp * Para.CDRAG * pi * a * sqrt(b^2 + c^2) *(d_bottom);       % Drag part of the work
% Work component for turbulence based on the shear velocity and the velocity gradient
u_star = 0.4 * sqrt(DynColl_Output(8)^2 + DynColl_Output(9)^2) / log(a/(2*0.001));  % mit z0=0.001 based on Cheng 1999
dvdz = 0.7 * sqrt(DynColl_Output(8)^2 + DynColl_Output(9)^2) / a;
W_t = 2/3 * pi * a * b * c * rho_bulk * u_star^2 * dvdz;                            % Turbulence part of the work
% Sum of work components
W_sum = W_f + W_d + W_t;
% Assigning the sum to the last entry of the output variable
DynColl_Output(10+k+1) = W_sum;                                                

% Entrainment calculation
surface = 2 * pi * a * sqrt(b^2+c^2);                                               % Surface area as seen in the vertical plane
V_entr = surface * Para.ALPHAB * sqrt(DynColl_Output(8)^2 + DynColl_Output(9)^2);   % Volume entrained per time step
dv = V_entr;                                                                        % Extra variable for accounting the overall volume change

%% Assign values to the output variable 
DynColl_Output(1) = amb(1);                                                         % 'u' of the cloud for calculating the x-position of the cloud centroid
DynColl_Output(2) = amb(2);                                                         % 'v' of the cloud for calculating the y-position of the cloud centroid
DynColl_Output(3) = - 3 * (a - abs(DynColl_Output(7)))/8;                           % 'w' of the cloud for calculating the z-position of the cloud centroid above ground 
DynColl_Output(4) = V_entr * rho_amb_tmp;                                           %  'dm/dt' rate of ambient fluid mass entrained 
% DynColl_Output(5) and (6) are the kinetic and potential energy that are changed in "main" for bookkeeping
% DynColl_Output(7), (8) and (9) have already been accounted earlier

%Equation for the bottom shear-stress from current only, including the Chézy Coefficient (here set to 65) at the end that can be taken from the far-field model
tau_b = (9.81 * rho_amb_tmp * amb(1)^2 + amb(2)^2) / 65^2;                          % A value of 65 is chosen for a sandy bed.

% Total sediment volume
Sed_mass = sum(Cloud_Input(11:10+k).*Sediment.rho');
% Homogeneous sediment concentration
C_solid = Sed_mass/volumen;
% Iteration over the sediment components for settling/deposition after Krone (1962)
for ik=1:k
    
    % Calculate settling velocity of cohesive sediment based on sediment concentration to check for hindered settling
    % (conservative approach: take the max. settling velocity to consider high settling velocities specified by the user)
    if Sediment.cohe(ik) == 1
        if C_solid <= 25
            Sediment.ws(ik) = max([0.000034*0.3048 Sediment.ws(ik)]);
        elseif C_solid > 3000
            Sediment.ws(ik) = max([0.0069 * 0.3048 Sediment.ws(ik)]);
        else
            Sediment.ws(ik) = max([(0.0000225+1.6*10^(-7)* C_solid^(4/3))*0.3048 Sediment.ws(ik)]);
        end
    end
    ws = abs(Sediment.ws(ik));                                                      % absolute settling velocity for later easier use
    cb = ws * dt * pi * b * c * Cloud_Input(10+ik) / volumen;                       % near bed sediment concentration for each fraction
    % Determine deposition flux reduction factor based on user input
    if Para.DELTA == -1
        f_reduction = 1 - tau_b / (Sediment.shear(ik) * 9.81);
        if f_reduction < 0
            f_reduction = 0.0;
        end
    else
        f_reduction = Para.DELTA;
    end
    % Calculate deposition flux
    V_setl = ws * cb * f_reduction;
    % Accounting the mass and volume reduction due to settling for each sediment fraction
    DynColl_Output(4) = DynColl_Output(4) - V_setl * Sediment.rho(ik);      	    % 'dm/dt' Step:2 rate of mass loss due to settling (Eq:3.1 Brandsma 1976)
    DynColl_Output(10+ik) = - V_setl;                                               % 'dP/dt' Solids that are lost per sediment fraction (Eq: 3.4 Brandsma 1976)
    dv = dv - V_setl;                                                               % 'dv' as change in total volume for bookkeeping
end
% Total change of the volume
DynColl_Output(10) = dv;                                                            

%% Return
if nargout < 1
    'not enough output arguments included when calling the function. Make sure that the cloud is returned'
elseif nargout == 1
    varargout{1} = DynColl_Output;
end
end

