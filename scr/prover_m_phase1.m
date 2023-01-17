function [varargout] = prover_m_phase1(Cloud_Input,Hopper,Sediment,Coef,d,Rho_amb,V_amb,z,iz,varargin)

%% prover_m_phase1: Convective Descent
% Made by:          Jannek Gundlach
% Date of Change:   13.01.2023
% Definition:       Sediment entering the water column and descending
%                   Calculating the change in mass, momentum, buoyancy, 
%                   vorticity and solids per Runge-Kutta 4 (RK4) fraction of the time step 
% Phase 1 
% End:              Bottom impact (checked at the end)

%% Define input   
% practical values if non defined
defaultg = 9.81;                % gravitational constant

% input parser
p = inputParser;
p.FunctionName = 'prover_m_phase1';
p.CaseSensitive = true;
p.StructExpand = true;
validStruct = @(x) isstruct(x);
validDouble = @(x) isnumeric(x) && isfloat(x);
validPosDouble = @(x) isnumeric(x) && isfloat(x) && all(x>=0);
validMatrix = @(x) ismatrix(x) && all(x>=0);
addRequired(p,'Cloud_Input',validDouble);
addRequired(p,'Hopper',validStruct);
addRequired(p,'Sediment',validStruct);
addRequired(p,'Coef',validStruct);
addRequired(p,'d',validPosDouble);
addRequired(p,'Rho_amb',validPosDouble);
addRequired(p,'V_amb',validDouble);
addRequired(p,'z',validMatrix);
addRequired(p,'iz',validPosDouble);
addOptional(p,'g',defaultg,validPosDouble);

parse(p,Cloud_Input,Hopper,Sediment,Coef,d,Rho_amb,V_amb,z,iz,varargin{:});

% Number of sediment fractions
k = length(Sediment.name(:));


%% Input preparation
% Check for valid input
if Cloud_Input(3) <= 0.0
    "Sediment cloud is not descending, please check your input and the definition of your z coordinate (Depth is defined positively)"
end
% iterate for the right depth coordinate and defined ambient density
while Cloud_Input(3) >= z(iz)
    iz = iz + 1;
end
% Interpolate ambient density at cloud depth 
Rho_amb_tmp = Rho_amb(iz) - (Cloud_Input(3) - z(iz)) * (Rho_amb(iz) - Rho_amb(iz-1)) / (z(iz) - z(iz-1));
% Density gradient equation 3.16 in Brandsma & Divoky 1976
drhodz = (Rho_amb(iz) - Rho_amb(iz-1)) / (z(iz) - z(iz-1)); 
% Calculate volume and cloud density
volumen = (Cloud_Input(4) + Cloud_Input(8)) /Rho_amb_tmp;
rho_bulk = Cloud_Input(4) / volumen;

% Check for correct ratio between ambient and cloud density
if rho_bulk <= Rho_amb_tmp
    "Sediment cloud will not reach the ground. The model is designed for the bottom impact. Some processes will be missing."
end
% Radius of the cloud
b = (1.5 * volumen / pi)^(1 / 3);




%% Calculate entrainment coefficient
alpha0 = Coef.ALPHA0;
alpha =  alpha0 * sqrt(tanh((p.Results.g * volumen * (rho_bulk - Rho_amb_tmp) / (2.0 * pi * 0.16 * rho_bulk * Cloud_Input(9)^2 * Coef.ALPHA0))^2));


%% Main calculations of phase 1

% Mass of the cloud including the mass coefficient
cm_mass = Coef.CM * Cloud_Input(4);
% Motion of the cloud in x-, y-, z-direction
uu = Cloud_Input(5) / cm_mass;
vv = Cloud_Input(6) / cm_mass;
ww = Cloud_Input(7) / cm_mass;
% Total motion of the cloud
phi = sqrt((uu - V_amb(1))^2 + (vv - V_amb(2))^2 + (ww)^2);

% Entrained volume
V_entr = 2.0 * pi * b^2 * alpha * phi;                      % (equation 3.5 in Brandsma 1976)
% Bookkeeping of the volume changes (entrainment adds)
dv = V_entr;

% Calculate the output of the cloud per time step interval based on RK4
ConvDesc_Output(1) = uu;                                    % 'dx/dt' Change in x-position of the cloud
ConvDesc_Output(2) = vv;                                    % 'dy/dt' Change in y-position of the cloud
ConvDesc_Output(3) = ww;                                    % 'dz/dt' Change in z-position of the cloud
ConvDesc_Output(4) = V_entr * Rho_amb_tmp;                  % 'dm/dt' Step:1 rate of ambient fluid mass entrained  
drag = Coef.CD * Rho_amb_tmp * pi * b^2 * phi * 0.5;        % Drag force (equation 3.2 in Brandsma 1976)
ConvDesc_Output(5) = V_entr * Rho_amb_tmp * V_amb(1) - drag * (uu - V_amb(1)) * 0.5; % 'dIx/dt' Step:1 Change of momentum in x- direction due to entrainment and drag
ConvDesc_Output(6) = V_entr * Rho_amb_tmp * V_amb(2) - drag * (vv - V_amb(2)) * 0.5; % 'dIy/dt' Step:1 Change of momentum in y- direction due to entrainment and drag
ConvDesc_Output(7) = volumen * (rho_bulk - Rho_amb_tmp) * p.Results.g - drag * ww;         % 'dIz/dt' Step:1 Change of momentum in z- direction due to entrainment and drag
ConvDesc_Output(8) = V_entr * (Rho_amb(1) - Rho_amb_tmp);                        % 'dB/dt' Step:1 Change of buoyancy due to entrainment
ConvDesc_Output(9) = -3.0 * b^2 * p.Results.g * drhodz / Rho_amb(1); % 'dW/dt' Change of vorticity due to density gradients (equation 3.14 in Brandsma 1976)
                                                            % 3.0 = vorticity-dissipation-coefficient according to Turner 1960

                                                        
%% Calculate sediment stripping for each sediment fraction

Stripping = 2.0 * pi * b^2 * phi * Coef.STRIP;              % Total amount of stripped material per time step fraction RK4

% Total sediment volume
Sed_mass = sum(Cloud_Input(11:end).*Sediment.rho');
% Homogeneous sediment concentration
C_solid = Sed_mass/volumen;

% Iterate stripping and settling of sediments for each sediment fraction
for ik=1:k
    
    % Calculate settling velocity of cohesive sediment based on sediment concentration to check for hindered settling
    % (conservative approach: take the max. settling velocity to consider high settling velocities specified by the user)
        if Sediment.cohe(ik) == 1
            % Calculate settling velocity based on sediment concentration
            if C_solid <= 25
                Sediment.ws(ik) = max([0.000034*0.3048 Sediment.ws(ik)]);
            elseif C_solid > 3000
                Sediment.ws(ik) = max([0.0069 * 0.3048 Sediment.ws(ik)]);
            else
                Sediment.ws(ik) = max([(0.0000225 + 1.6 * 10^(-7) * C_solid^(4/3))*0.3048 Sediment.ws(ik)]);
            end
        end
    
    ws = abs(Sediment.ws(ik));                          % defined settling velovity
    if ws <= ww
        beta_neu = 1.0;                                 % if the settling velocity is less than the downward cloud motion, no settling occurs
    else
        beta_neu = Coef.BETA;                           % else calculated according to the coefficient beta
    end
    % Calculate the amount of stripping for each sediment fraction that is able to be stipped
    Strip(ik) = Stripping * Cloud_Input(10+ik) / sum( Sediment.strip .* Cloud_Input(11:10+k)' );
    % Calculate the amount of settling for each sediment fraction
    V_setl = pi * b^2 * abs(Sediment.ws(ik)) * (1.0 - beta_neu) * Cloud_Input(10+ik) / volumen; 
    ConvDesc_Output(4) = ConvDesc_Output(4) - (V_setl+Strip(ik)) * Sediment.rho(ik);      % 'dm/dt' Step:2 Change of mass due to stripping and settling
    ConvDesc_Output(5) = ConvDesc_Output(5) - (V_setl+Strip(ik)) * Sediment.rho(ik) * uu; % 'dIx/dt' Step:2 Change of Momentum due to stripping and settling (Eq:3.2 Brandsma 1976)
    ConvDesc_Output(6) = ConvDesc_Output(6) - (V_setl+Strip(ik)) * Sediment.rho(ik) * vv; % 'dIy/dt' Step:2 Change of Momentum due to stripping and settling (Eq:3.2 Brandsma 1976)
    ConvDesc_Output(7) = ConvDesc_Output(7) - (V_setl+Strip(ik)) * Sediment.rho(ik) * ww; % 'dIz/dt' Step:2 Change of Momentum due to stripping and settling (Eq:3.2 Brandsma 1976)
    ConvDesc_Output(8) = ConvDesc_Output(8) - (V_setl+Strip(ik)) * (Rho_amb(1) - Sediment.rho(ik)); % EP(8) 'dB/dt' Step:2 Change in buoyancy due to stripping and settling (Eq:3.3 Brandsma 1976)
    % Bookkeeping of volume changes
    dv = dv - (V_setl+Strip(ik));                       % Stripping and Settling are accounted negatively
    ConvDesc_Output(10+ik) = -(V_setl+Strip(ik));       % 'dP/dt' Change in the amount of solids for each sediment fraction (Eq: 3.4 Brandsma 1976)
end
% Total change of the volume
ConvDesc_Output(10) = dv;

%% Return

if nargout < 1
    'not enough output arguments included when calling the function. Make sure that the cloud output is returned'
elseif nargout == 1
    varargout{1} = ConvDesc_Output;
elseif nargout >= 2
    varargout{1} = ConvDesc_Output;
    varargout{2} = iz;
end
end