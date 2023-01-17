%% Das ist der Diskretisierungssolver für Phase 1 & 2 nach dem Klassischen Runge Kutta Schema (4-Schrittig)
function [varargout] = prover_m_rk4(Cloud,b_impact,Hopper,Sediment,Para,d,dt,qe,amb,istep,z,iz,Rho_amb,varargin)
%% prover_m_phase1: Convective Descent
% Made by:          Jannek Gundlach
% Date of Change:   13.01.2023
% Definition:       Numerical solver based on the classic Runge-Kutta 4 (RK4) scheme.
%                   It calls the Convective Descent function and the Dynamic Collapse,
%                   dependent on the vertical position of the cloud. RK4 splits the numerical
%                   time step into sections that are weighted at the end and summed to get a
%                   more accurate result. For its definition check Kutta (1901) or check
%                   the Mathwork Course on "Solving ODEs in MATLAB", 3: Classical Runge-Kutta, ODE4. 

    defaultRhoAmb = [1030 1030];
    defaultAmb = [0.0 0.0];
    defaultiz = 1;
    defaultz = [0 30];
    defaultfbed = 0.0;
    defaulta0 = 10;
    
    p = inputParser;
        p.FunctionName = 'prover_m_rk4';
        p.CaseSensitive = true;
        p.StructExpand = true;
    validStruct = @(x) isstruct(x);
    validDouble = @(x) isnumeric(x) && isfloat(x);
    validPosDouble = @(x) isnumeric(x) && isfloat(x) && all(x>=0);
    validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
    validScalarNum = @(x) isnumeric(x) && isscalar(x);
    validMatrix = @(x) ismatrix(x) && all(x>=0);
    addRequired(p,'Cloud',validDouble);
    addRequired(p,'b_impact',validScalarNum);
    addRequired(p,'Hopper',validStruct);
    addRequired(p,'Sediment',validStruct);
    addRequired(p,'Para',validStruct);
    addRequired(p,'d',validScalarPosNum);
    addRequired(p,'dt',validDouble);
    addRequired(p,'qe',validDouble);
    addRequired(p,'amb',validDouble);
    addRequired(p,'istep',validScalarPosNum);
    addRequired(p,'z',validMatrix);
    addOptional(p,'iz',defaultiz,validPosDouble);
    addOptional(p,'Rho_amb',defaultRhoAmb,validMatrix);
    addOptional(p,'f_bed',defaultfbed,validPosDouble);
    addOptional(p,'a0',defaulta0,validScalarNum);
    
    parse(p,Cloud,b_impact,Hopper,Sediment,Para,d,dt,qe,amb,istep,z,iz,varargin{:});
    %% Automatic choice, if Convective Descent or Dynamic Collapse has to be called
    if b_impact == 1
        funktion = 'phase2';
    elseif b_impact == 0
        funktion = 'phase1';
    else
        error('Error! Cloud is neither in the first nor in the second phase. Check input and parameter b_impact')
    end
    
    %% Input preperation for Runge-Kutta 4th order
    rk4_inp = {Cloud zeros(size(Cloud)) zeros(size(Cloud)) zeros(size(Cloud)) zeros(size(Cloud))};
    rk4_fakt = {[0.5 0 0 0] [0 0.5 0 0] [0 0 1 0] [1 2 2 1]};
    rk4_var = zeros(4,length(Cloud));
    iz = iz;

    %% Execution of the four-step discretization according to classic Runge-Kutta
    for rk4 = 1:4
        % Switch between cases to call the right function 
        switch funktion
            case 'phase1'
                [ph_out,iz] = prover_m_phase1(rk4_inp{rk4},Hopper,Sediment,Para,d,p.Results.Rho_amb,amb,z,iz);
            case 'phase2'
                [ph_out] = prover_m_phase2(rk4_inp{rk4},Hopper,Sediment,Para,d,p.Results.Rho_amb,istep,qe,amb,z,dt,iz);
        end
        % Iteration over the equations for each time step in pahse 1 and 2 respectively
        for equ = 1:length(ph_out)  
            % Bookkeeping of the changes per RK4-step (rk4) and each cloud variable (equ) fr each time step
            rk4_var(rk4,equ) = dt * ph_out(equ);                            
            % If condition for getting the correct quotient in eacht RK4-step.
            if sum(rk4_fakt{rk4})>1
                quotient = sum(rk4_fakt{rk4});
            else
                quotient = 1;
            end
            rk4_inp{rk4+1}(equ) = rk4_inp{1}(equ) + sum(rk4_var(:,equ)' .* rk4_fakt{rk4}) / quotient; % combine the results to the next RK4-step
        end
    end
    varargout{1} = rk4_inp{5};
    varargout{2} = iz;