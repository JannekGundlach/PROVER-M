function [E1] = prover_m_main(app,file_info,d,amb,rho_amb_init,z_i,u_amb_crit,mass,draft_predisp,draft_postdisp,n_dump,Para,Sediment,dt_sim,t_sim_istep)
%% prover_m_phase1: Convective Descent
% Made by:          Jannek Gundlach and Maximilian Behnke
% Date of Change:   17.01.2023
% Definition:       Main script that is called by the GUI.
%                   Here, bookkeeping of the sediment cloud, sediment
%                   losses due to stripping and settling is done.
%                   During the Convective Descent the function of the first phase is
%                   called through the fourth order Runge Kutta scheme (RK4).
%                   During the Dynamic Collapse the function of the second
%                   phase is called through RK4.
%
% End:              The model terminates when the sediment cloud becomes
%                   passive, e.g. the ambient currents will have a stronger
%                   influence on the clouds motion than the spreading and motion
%                   due to the disposal.

% clean up environment
clearvars -except app file_info d amb rho_amb_init z_i u_amb_crit mass draft_predisp draft_postdisp n_dump Para Sediment dt_sim t_sim_istep
close all
clc

% Read & assign variables
close all                                                                       % close all previous figures
textfilename = app.InputcaseDropDown.Value;                                     % name of simulation
name = char(strcat(textfilename(1:(end-4))));                                   % pure name of control file for later use in saving output

Hopper.draft_predisp = draft_predisp;                                           % Pre-disposal Draft of the Hopper [m]
Hopper.draft_postdisp = draft_postdisp;                                         % Post-disposal Draft of the Hopper [m]

mass_dump = mass/n_dump;                                                        % Total mass proportion per dumping instance
n_rho_amb = length(rho_amb_init);                                               % # of vertical layers
k = length(Sediment.name);                                                      % # of solid fractions
g = 9.81;                                                                       % Gravitational acceleration [m/s]
diskretisierung = 'CRK';                                                        % Discretization scheme (CRK - Classical Runge Kutta)

for i_dump = 1:n_dump                                                           % loop over individual dumping instances

    clear Cloud VLOSS E_counter t                                               % clear reoccuring output variables

    fid = fopen(strcat(file_info{1},'/output/',name,'_dump',num2str(i_dump),'.txt'),'w');             % open output file for writing cloud results
    fid2 = fopen(strcat(file_info{1}, '/output/',name,'_VLOSS_dump',num2str(i_dump),'.txt'),'w');     % open output file for writing lossed material results (VLOSS)

    % preallocation
    [drho, VLOSS, E_pot, E_kin, KE_P2, PE_P2, delta_Work, delta_KE, delta_PE] = deal([]);   % for growing vectors in every loop
    [b, u, v, w, t, x, y, z, m, V, rho, W, a, c, sed, Sed, fC, Work, KE, PE] = deal([]);    % for the Cloud structure
    Cloud = table(b, u, v, w, t, x, y, z, m, V, rho, W, a, c, sed, Sed, fC, Work, KE, PE);  % create a table with the matching Cloud.i variables
    Cloud = table2struct(Cloud); % convert to structure

    Cloud(1).b = ( mass_dump / (2.0 * pi / 3.0))^(1 / 3);                       % Initial radius of the cloud

    %% Start of Disposal: Initiation of Phase 1 (Convective Descent)
    % give feedback to command window and in the App
    feedText = {'Initial conditions are specified'};                            % feedback in the GUI command feed
    app.commandfeedTextArea.Value = feedText;
    disp('"Initial conditions are specified"')                                  % feedback in the command window


    %% Calculation
    RB = Cloud(1).b;                                                            % save initial cloud radius
    rho_bulk = (sum(Sediment.rho .* Sediment.frac) + (1 - sum(Sediment.frac))...% bulk density of sediments-water mixture
        * mean(rho_amb_init));

    % Initialization                                                            % Abbrevation | Description
    Cloud(1).u = 0.0;                                                           %     'u'     | initial cloud velocity in x-direction
    Cloud(1).v = 0.0;                                                           %     'v'     | initial cloud velocity in y-direction
    Cloud(1).w = sqrt((2*8*g - 2*(1/2.65)*g*5)/(1+0.02*g*rho_bulk));            %     'w'     | initial cloud velocity in z-direction


    % Definition and calculation of the ambient density
    rho_amb = rho_amb_init;                                                     % Definition of the momentary (near-surface) density


    %Definition of the initial conditions of the sediment cloud                 Abbreviation   | Description
    Cloud_Input(1) = 0.0;                                                       %     'x'       | Initial position of the cloud in x-direction from disposal location
    Cloud_Input(2) = 0.0;                                                       %     'y'       | Initial position of the cloud in y-direction from disposal location
    Cloud_Input(3) = Hopper.draft_postdisp + ...
        (Hopper.draft_predisp - Hopper.draft_postdisp)/2 + 3/8 * RB ;           %     'z'       | Initial position of the cloud in z-direction from disposal location
    volumen = (4.0 / (2.0 * 3.0)) * pi * Cloud(1).b^3;                          %     'V'       | Initial volume of the cloud
    Cloud_Input(4) = rho_bulk * volumen;                                        %     'm'       | Initial mass of the cloud
    cm_mass = Para.CM * Cloud_Input(4);                                         %     'cm'      | Initial mass portion of the total impulse; 1,0 < CM < 1,5 (Brandsma 1976)
    Cloud_Input(5) = cm_mass * Cloud(1).u;                                      %     'Ix'      | Initial impulse in x-direction
    Cloud_Input(6) = cm_mass * Cloud(1).v;                                      %     'Iy'      | Initial impulse in y-direction
    Cloud_Input(7) = cm_mass * Cloud(1).w;                                      %     'Iz'      | Initial impulse in z-direction
    Cloud_Input(8) = (rho_amb_init(1) - rho_bulk) * volumen;                    %     'B'       | Initial (negative) buoyancy due to the density difference
    Cloud_Input(9) = Cloud(1).b * Cloud(1).w * Para.KV;                         %     'W'       | Initial vorticity

    Cloud_Input(10) = volumen;                                                  %  Initial volume of the cloud

    fluid_V = volumen;                                                          %     V_f       | Temporary variable for the determination of the fluid volume (w/out solid fractions)
    for i=1:k                                                                   % loop variable i over all solid fractions
        Cloud_Input(10+i) = Sediment.frac(i) * volumen;                         %     Sed_V     | store the absolute fraction volumes
        fluid_V = fluid_V - Cloud_Input(10+i);                                  % Bookkeeping of the remaining volume
    end

    clear i;                                                                    % delete loop variable to avoid doubling
    t(1) = 0.0;                                                                 % set initial time
    iz = 1;                                                                     % iteration variable for the vertical layer
    istep = 1;                                                                  % set initial # of time step
    b_impact = 0;                                                               % bottom encountered? [0/1] if 1 dynamic collapse starts

    % Time step definition dt
    Delta0 = (rho_bulk-rho_amb_init(1)) / rho_amb_init(1);                      % Initial (negative) buoyancy
    E1 = (rho_amb_init(n_rho_amb) - rho_amb_init(1)) / (d * rho_amb_init(1));   % Vertical, relative density difference
    FF = Cloud(1).w / sqrt(g * Delta0 * RB);                                    % Dimensionless descent rate/velocity (w~ after Escudier & Maxworthy 1973 w/out entrainment coefficient)

    dt_uni = 0.01*Cloud(1).b*FF*((1.0+Para.ALPHA0*d/Cloud(1).b)^2)...           % time step definition for a well-mixed water column (after Koh & Chang 1973)
        /(Cloud(1).w*2.0);
    dt_strat = 0.01*pi*FF/(abs(E1 * RB / Delta0)^0.5);                          % time step definition for a stratified water column (after Koh & Chang 1973)
    if Cloud(1).w == 0.0                                                        % check conditions to decide on the time step definition
        dt = 1.0;
    else
        if rho_amb_init(end) ~= rho_amb_init(1)                                 % corresponds to stratification
            dt = min(dt_strat,dt_uni);
        else
            dt = dt_uni;                                                        % corresponds to well-mixed water column
        end
        % adjust time step to desired far-field time step if in reasonable deviation
        frac_dt = dt_sim/round(dt_sim/dt);                                      % find nearest fraction of dt_sim
        dev_dt = max(dt,frac_dt)/min(dt,frac_dt) - 1;                           % deviation of dt from the nearest fraction
        if dev_dt <= 0.2                                                        % only adjust the time step if the deviation is lower than or equal to 20%
            dt = frac_dt;
        else
            dt = dt;
            feedText = [feedText; 'Desired output time step interval could not be realized. Consider adjusting the value.']; 
            app.commandfeedTextArea.Value = feedText;                           % give feedback to command feed
        end
    end

    %% Time iteration in Phase 1 (Convective Descent)
    % Give feedback in the command window
    feedText = [feedText; 'Start calculation/time iteration'; 'Begin with the calculation of Phase 1'];
    app.commandfeedTextArea.Value = feedText;                                   % feedback in the GUI command feed
    drawnow                                                                     % update the text in the live feed

    disp('"Start calculation/time iteration"')                                  % feedback in the command window
    disp('"Begin with the calculation of Phase 1"')                             % feedback in the command window
    % Write headers in the textfile for cloud results
    fprintf(fid, ['Disposal process ', num2str(i_dump), ' of ', num2str(n_dump),'\n']);
    fprintf(fid, 'Begin with the calculation of Phase 1\n');
    format_tmp_Sed = cell(1,length(Sediment.name));
    format_tmp_Sed(:)={'%8.3f |'};
    formatSpec = strcat('%8.4f | %5.2f | %5.2f | %5.2f | %4.2f | %4.2f | %4.2f | %6.2f | %6.2f | %6.2f | ', strjoin(string(format_tmp_Sed)),' %11.2f | %8.2f | %8.2f | %4.2f | %6.4f \n'); % specify column properties for the output file
    fprintf(fid, strcat('Time [s] | x [m] | y [m] | z [m] |  u &  v &  w [m/s] |  a [m] |  b [m] |  c [m] | Sand[m³]| Silt[m³] | Clay[m³] | Masse [kg]  | Vol [m³] |rho[kg/m³]| Wirb | fluidC \n')); % write column headers in output file
    fprintf(fid, '-------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n'); % seperation line in output file

    % Write headers in the textfile for lossed material (VLOSS)
    formatSpecVL = strcat('%8.4f | %5.2f | %5.2f | %5.2f | %9.2f | %9.2f | %8.2f | %10.4f | %6.2f | %6.2f | \n'); % specify column properties for the output file
    fprintf(fid2, strcat('Time [s] | x [m] | y [m] | z [m] | Loss of Sand, Silt, Clay in [m³] | fluid in %% |  a [m] |  b [m] | \n')); % write column headers in output file
    fprintf(fid2, '---------------------------------------------------------------------------------------------------- \n'); % seperation line in output file
    fprintf(fid2, 'Material loss due to stripping \n'); % seperation line in output file
    fprintf(fid2, '---------------------------------------------------------------------------------------------------- \n'); % seperation line in output file

    t_final = 1800;                                                             % maximum simulation time definition [s]
    t_end_ph1 = t_final;                                                        % store final simulation time for phase 1
    cloud_counter_phase1 = 1;                                                   % for the export of turbidity clouds at UNTRIM time steps (dt_sim)
    idx_alt = 1;                                                                % for the export of turbidity clouds at UNTRIM time steps (dt_sim)

    % Calculation loop for Phase 1
    while t(istep) <= t_end_ph1                                                 % loop while below the maximum simulation time (maximum run time)

        Cloud(istep).t = t(istep);                                              % Description
        Cloud(istep).x = Cloud_Input(1);                                        % position of the cloud in x-direction from disposal location
        Cloud(istep).y = Cloud_Input(2);                                        % position of the cloud in y-direction from disposal location
        Cloud(istep).z = Cloud_Input(3);                                        % position of the cloud in z-direction (corresponds to bottom edge of hopper dredger)
        Cloud(istep).m = Cloud_Input(4);                                        % mass of the cloud
        Cloud(istep).V = Cloud_Input(10);                                       % volume of the cloud
        qe = 0.0;                                                               % rate of volume change
        volumen = (Cloud_Input(4) + Cloud_Input(8)) / rho_amb_init(1);          % Volume of the cloud
        cm_mass = 1.0 / (Para.CM * Cloud_Input(4));                             % mass portion of the total impulse; 1,0 < CM < 1,5 (Koh & Chang 1973; Brandsma & Divoky 1976)
        Cloud(istep).u = Cloud_Input(5) * cm_mass;                              % Velocity of the cloud in x-direction
        Cloud(istep).v = Cloud_Input(6) * cm_mass;                              % Velocity of the cloud in y-direction
        Cloud(istep).w = Cloud_Input(7) * cm_mass;                              % Velocity of the cloud in z-direction
        if Cloud_Input(9) <= 0.0                                                % reset vorticity if needed
            Cloud_Input(9) = 0.0;
        end
        Cloud(istep).rho = rho_bulk;                                            % bulk density of sediments-water mixture
        Cloud(istep).W = Cloud_Input(9);                                        % Vorticity
        Cloud(istep).a = (1.5 * volumen / pi) ^(1/3);                           % vertical radius of the cloud in [m]
        Cloud(istep).b = Cloud(istep).a;                                        % horizontal diameter of the cloud in [m]
        Cloud(istep).c = 0.0;                                                   % horizontal radius of the cloud in [m] (Dummy in phase 1)
        fluid_V = Cloud_Input(10);
        for ik=k:-1:1                                                           % loop variable i over all solid fractions
            sedK(ik) = Cloud_Input(10+ik) / volumen;                            % store the relative fraction volumes
            SedK(ik) = Cloud_Input(10+ik);                                      % store the absolute fraction volumes
            fluid_V = fluid_V - Cloud_Input(10+ik);                             % Bookkeeping of the remaining volume
        end
        Cloud(istep).sed = sedK;                                                % Assign initial relative fraction volumes
        Cloud(istep).Sed = SedK;                                                % Assign initial absolute fraction volumes
        Cloud(istep).fC = fluid_V / volumen;                                    % fluid volume of the cloud

        Cloud(istep).Work = 0;                                                  % Initial work components (for Phase 2 Dynamic Collapse)

        % Calculation of the bulk density of the sediment cloud
        rho_bulk = Cloud_Input(4) / volumen;

        % Interpolation of the ambient density at the cloud center of gravity
        rho_amb_tmp = rho_amb(iz) + (Cloud_Input(3) - z_i(iz)) * (rho_amb(2) - rho_amb(1)) / (z_i(2) - z_i(1));

        % Calculation of the density difference
        drho(istep) = rho_bulk - rho_amb_tmp;


        %% Check if cloud reached the bottom
        if (Cloud(istep).z + 5.0 * Cloud(istep).a / 8.0) >= d
            b_impact = 1;                                                       % bottom reached
        else
            b_impact = 0;                                                       % bottom not reached
            fprintf(fid, formatSpec, t(istep),Cloud(istep).x, Cloud(istep).y, Cloud(istep).z, Cloud(istep).u, Cloud(istep).v, Cloud(istep).w, Cloud(istep).a, Cloud(istep).b, Cloud(istep).c, Cloud(istep).Sed, Cloud(istep).m, Cloud(istep).V, Cloud(istep).rho, Cloud(istep).W, Cloud(istep).fC); % print current cloud results in output file
        end

        % Initialization of the bottom encounter (Begin Phase 2)
        if b_impact == 1                                                        % bottom reached
            if istep ==1                                                        % bottom reached instantly
                t_end_ph1 = t(istep);                                           % store final time step of phase 1
                fprintf(fid, 'Water depth too shallow for Phase 1 - direct transition to Phase 2\n');   % write in output file
            else
                t_end_ph1 = t(istep-1);                                         % store final time step of phase 1
                fprintf(fid, '-------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n'); % seperation line in output file
            end
            continue
        else

            %% Calculation of the dissolved/stripped material at every far field model-time step dt_sim
            if dt_sim * cloud_counter_phase1 <= t(istep)
                C1 = 2.0 * pi / 3.0;                                            % geometrical factor for later calculations
                % Determine output variables
                VLOSS(cloud_counter_phase1,1) = t(istep);                       % time
                VLOSS(cloud_counter_phase1,2) = Cloud(istep).x;                 % position of the cloud in x-direction
                VLOSS(cloud_counter_phase1,3) = Cloud(istep).y;                 % position of the cloud in y-direction
                VLOSS(cloud_counter_phase1,4) = Cloud(istep).z;                 % position of the cloud in z-direction
                VLOSS(cloud_counter_phase1,k+5) = Cloud(istep).fC;              % fluid volume of the cloud
                VLOSS(cloud_counter_phase1,k+6) = Cloud(istep).a;               % vertical radius of the cloud
                VLOSS(cloud_counter_phase1,k+7) = Cloud(istep).b;               % horizontal diameter of the cloud
                for k1 = 1:k
                    VLOSS(cloud_counter_phase1,k1+4) = Cloud(idx_alt).Sed(k1) - Cloud(istep).Sed(k1);
                end

                fprintf(fid2, formatSpecVL, VLOSS(cloud_counter_phase1,1), VLOSS(cloud_counter_phase1,2), VLOSS(cloud_counter_phase1,3), VLOSS(cloud_counter_phase1,4), VLOSS(cloud_counter_phase1,5), VLOSS(cloud_counter_phase1,6), VLOSS(cloud_counter_phase1,7), VLOSS(cloud_counter_phase1,8), VLOSS(cloud_counter_phase1,9), VLOSS(cloud_counter_phase1,10)); % print current vloss results in output file
                cloud_counter_phase1 = cloud_counter_phase1 + 1;                % increase the counter
                idx_alt = istep;                                                % store current time step
            end
        end

        %% Updating Live GUI parameters for Convective Descent
        if t_sim_istep == 1 && mod(dt*(istep-1),dt_sim) == 0 && t(istep) > 0 ...
                || t_sim_istep == 0 && mod(istep,25) == 0                       % update Live parameters in the GUI for Phase 1 for every xth (25) step

            % Set axis limits of live plot in GUI
            app.UIAxes.YLim = [-app.WaterdepthEditField.Value 0];
            app.UIAxes.YLabel.String = 'Water depth [m]';
            app.UIAxes.XLim = [0 app.WaterdepthEditField.Value]; 

            % update parameters
            app.TimeEditField.Value = round(Cloud(istep).t,1);                  % update current simulation time
            app.StrippingEditField.Value = sum(sum(VLOSS(1:cloud_counter_phase1-1,5:(k+4))));   % update stripping volume
            app.CloudradiusEditField.Value = round(Cloud(istep).a,1);           % update cloud radius
            drawnow                                                             % show updated parameters in GUI

            % plot cloud contour
            [x_sphere,y_sphere,z_sphere] = sphere;                              % create a sphere
            r_ph1_plot = Cloud(istep).b;                                        % assign radius of the sphere
            x_half_sphere = x_sphere(1:11,[1 11]).* r_ph1_plot;                 % get x-coordinates of the cloud sphere
            y_half_sphere = y_sphere(1:11,[1 11]).* r_ph1_plot;                 % get y-coordinates of the cloud sphere
            z_half_sphere = z_sphere(1:11,[1 11]).* r_ph1_plot;                 % get z-coordinates of the cloud sphere
            z_top_ph1_plot = Cloud(istep).z - 3/8 * r_ph1_plot;                 % get centre of mass coordinate

            hold(app.UIAxes,"off");                                             % axes setting
            plot(app.UIAxes, x_half_sphere,(z_half_sphere-z_top_ph1_plot),'Color','k')      % plot cloud contour
            drawnow
        end

        %% Calculation of the gradients with the numerical solver
        [Cloud_Input,iz] = prover_m_rk4(Cloud_Input,b_impact,Hopper,Sediment,Para,d,dt,qe,amb,istep,z_i,iz,rho_amb); % enter numerical solver (CRK)

        %% prepare next time step
        t(istep+1)=istep*dt;                                                    % time for next time step
        istep=istep+1;                                                          % increase the time step counter
        
    end

    %% Time iteration in Phase 2 (Dynamic Collapse)

    % Preparation of Phase 2 & transition in calculation scheme (equations of motion (Phase 1) to energy equations (Phase 2))
    fprintf(fid2, '---------------------------------------------------------------------------------------------------- \n'); % seperation line in output file
    fprintf(fid2, 'Material loss due to settling \n');                          % seperation line in output file
    fprintf(fid2, '---------------------------------------------------------------------------------------------------- \n'); % seperation line in output file

    p1_step_p2 = istep;                                                         % assign time step for transition between Convective Descent and Dynamic Collapse

    v_ges = Cloud(istep).u^2 + Cloud(istep).v^2 + Cloud(istep).w^2;             % Cloud velocity (for kinetic energy calculation)
    E_kin_P1 = (1/3) * rho_bulk * pi * Cloud(istep).a^3 * abs(v_ges);           % Initial kinetic energy
    E_pot_P1 = (3/8) * (Cloud(istep).rho - rho_amb_tmp) * g * Cloud(istep).V * Cloud(istep).a;   % Initial potential energy

    Cloud(istep).a = Cloud(istep).a;                                            % assign vertical cloud dimension
    Cloud(istep).c = Cloud(istep).b;                                            % assign horizontal cloud dimension
    Cloud(istep).b = Cloud(istep).b;                                            % assign horizontal cloud dimension
                                                                                % Abbreviation | Description
    Cloud_Input(1) = Cloud(istep).x;                                            %   'x'        | Initial position of the cloud in x-direction from disposal location
    Cloud_Input(2) = Cloud(istep).y;                                            %   'y'        | Initial position of the cloud in y-direction from disposal location
    Cloud_Input(3) = Cloud(istep).z;                                            %   'z'        | Initial position of the cloud in z-direction (corresponds to bottom edge of hopper dredger)
    Cloud_Input(5) = E_kin_P1;                                                  % 'Ekin'       | Initial kinetic energy
    Cloud_Input(6) = E_pot_P1;                                                  % 'Epot'       | Initial potential energy
    Cloud_Input(7) = Cloud(istep).a;                                            %   'a'        | Initial dimension of the cloud in z-direction
    Cloud_Input(8) = Cloud(istep).b;                                            %   'b'        | Initial dimension of the cloud in x-direction
    Cloud_Input(9) = Cloud(istep).c;                                            %   'c'        | Initial dimension of the cloud in y-direction
    Cloud_Input(11+k) = 0;                                                      % 'Work'       | Initial work component for the kinetic energy
    Cloud(istep).KE = E_kin_P1;                                                 % assign kinetic energy as a cloud variable
    Cloud(istep).PE = E_pot_P1;                                                 % assign potential energy as a cloud variable

    % Write results of the transitional ellipsoid in output file
    fprintf(fid, ' Ellipsoid for the transition from Phase 1 to Phase 2 \n');   % transitional ellipsoid
    fprintf(fid, ' -------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n'); % seperation line in output file
    fprintf(fid, formatSpec, t(istep),Cloud(istep).x, Cloud(istep).y, Cloud(istep).z, Cloud(istep).u, Cloud(istep).v, Cloud(istep).w, Cloud(istep).a, Cloud(istep).b, Cloud(istep).c, Cloud(istep).Sed, Cloud(istep).m, Cloud(istep).V, Cloud(istep).rho, Cloud(istep).W, Cloud(istep).fC); % print current cloud results in output file
    fprintf(fid, ' -------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n'); % seperation line in output file
    %% Begin Phase 2 - Calculation
    % Give feedback in the command window
    feedText = [feedText; 'Begin with the calculation of Phase 2'];
    app.commandfeedTextArea.Value = feedText;                                   % feedback in the GUI command feed
    drawnow                                                                     % update the text in the live feed

    disp('"Begin with the calculation of Phase 2"')                             % feedback in the command window
    fprintf(fid, ' Begin with the calculation of Phase 2 \n');                  % feedback in the command window
    fprintf(fid, ' -------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n'); % seperation line in output file

    t_end = t_final;                                                            % maximum simulation time definition (same as for phase 1)
    cloud_counter_phase2 = cloud_counter_phase1;                                % for the export of turbidity clouds at far-field model time steps (dt_sim)

    % Calculation loop for Phase 2
    while t(istep) <= t_end                                                     % loop while below the maximum simulation time (maximum run time)

        if p1_step_p2 == istep                                                  % Query to avoid double allocation for the transition time step
            E_counter = 1;                                                      % counter to keep track of time steps in Phase 2
            PE_P2(E_counter) = E_pot_P1;                                        % assign potential energy (2nd attempt)
            KE_P2(E_counter) = E_kin_P1;                                        % assign kinetic energy (2nd attempt)

            E_counter = E_counter + 1;                                          % increase the counter

            % prepare next time step
            t(istep+1)=istep*dt;                                                % time for next time step
            istep=istep+1;                                                      % increase the time step counter

            continue

        else
            %% Update Cloud-variables from CRK-Cloud_Input-calculation
            Cloud(istep).t = t(istep);                                          % time
            Cloud(istep).x = Cloud_Input(1);                                    % position of the cloud in x-direction from disposal location
            Cloud(istep).y = Cloud_Input(2);                                    % position of the cloud in y-direction from disposal location
            Cloud(istep).z = d - (3 * Cloud_Input(7) / 8);                      % position of the cloud in z-direction (corresponds to bottom edge of hopper dredger)
            Cloud(istep).m = Cloud_Input(4);                                    % mass of cloud
            Cloud(istep).V = Cloud_Input(10);                                   % volume of cloud
            Cloud(istep).rho = Cloud(istep).m / Cloud(istep).V;                 % density of cloud
            qe = Cloud(istep).V - Cloud(istep-1).V;                             % rate of volume change due to entrainment
            Cloud(istep).u = Cloud(istep-1).u;                                  % velocity of the cloud in x-direction
            Cloud(istep).v = Cloud(istep-1).v;                                  % velocity of the cloud in y-direction
            if istep == p1_step_p2 + 1
                Cloud(istep).w = Cloud(istep-1).w;                              % velocity of the cloud in z-direction
            else
                Cloud(istep).w = abs((Cloud(istep).z -Cloud(istep-1).z)) / dt;  % velocity of the cloud in z-direction
            end
            Cloud(istep).a = Cloud_Input(7);                                    % dimension of the cloud in z-direction
            Cloud(istep).b = Cloud_Input(8);                                    % dimension of the cloud in x-direction
            Cloud(istep).c = Cloud_Input(9);                                    % dimension of the cloud in y-direction
            Cloud(istep).Work = Cloud_Input(11+k);                              % work component for the kinetic energy
            Cloud_Input(11+k) = 0;                                              % reset work components to 0 for next iteration
            % Determine kinetic and potential energy
            PE_P2(E_counter) = (3/8) * (Cloud(istep).rho - rho_amb_tmp) * g * Cloud(istep).V * Cloud(istep).a;      % potential energy
            delta_PE(E_counter) = PE_P2(E_counter) - PE_P2(E_counter-1);        % difference in potential energy
            delta_Work(E_counter) = -abs(Cloud(istep).Work);                    % calculate difference in work
            if delta_Work(E_counter) == 0                                       % case when no work is done
                delta_KE(E_counter) = 0.0;                                      % no change in kinetic energy
            else                                                                % case when work is done
                delta_KE(E_counter) = delta_Work(E_counter) - delta_PE(E_counter);      % get change in kinetic energy from potential energy and work
            end
            mass_factor = (Cloud(istep).m/Cloud(istep-1).m);                    % accounting of the effect of entrainment and settling on the kinetic energy
            KE_P2(E_counter) = (KE_P2(E_counter-1) + delta_KE(E_counter))/mass_factor;        % kinetic energy at the current (E-counter) time step
            Cloud_Input(5) = KE_P2(E_counter);                                  % assign kinetic energy as Input for the next time step
            Cloud(istep).KE = KE_P2(E_counter);                                 % assign kinetic energy as a cloud variable
            Cloud(istep).PE = PE_P2(E_counter);                                 % assign potential energy as a cloud variable

            E_counter = E_counter + 1;                                          % increase counter
            fluid_V = Cloud(istep).V;                                           % assign fluid volume

            for ik=1:k                                                          % Loop variable i over all solid fractions
                sedK(ik) = Cloud_Input(10+ik) / Cloud(istep).V;                 % store the relative fraction volumes
                SedK(ik) = Cloud_Input(10+ik);                                  % store the absolute fraction volumes
                fluid_V = fluid_V - Cloud_Input(10+ik);                         % Bookkeeping of the remaining volume
            end
            % assign determined volume variables to the cloud variables
            Cloud(istep).sed = sedK;                                            % relative fraction volumes
            Cloud(istep).Sed = SedK;                                            % absolute fraction volumes
            Cloud(istep).fC = fluid_V / Cloud(istep).V;                         % fluid volume of the cloud
            fprintf(fid, formatSpec, t(istep),Cloud(istep).x, Cloud(istep).y, Cloud(istep).z, Cloud(istep).u, Cloud(istep).v, Cloud(istep).w, Cloud(istep).a, Cloud(istep).b, Cloud(istep).c, Cloud(istep).Sed, Cloud(istep).m, Cloud(istep).V, Cloud(istep).rho, Cloud(istep).W, Cloud(istep).fC); % print current cloud results in output file

        end

        %% Calculation of the settled material at every far-field time step dt_sim
        if dt_sim * cloud_counter_phase2 <= t(istep)
            C1 = 2.0 * pi / 3.0;
            % Determine output variables
            VLOSS(cloud_counter_phase2,1) = t(istep);
            VLOSS(cloud_counter_phase2,2) = Cloud(istep).x;
            VLOSS(cloud_counter_phase2,3) = Cloud(istep).y;
            VLOSS(cloud_counter_phase2,4) = Cloud(istep).z;
            VLOSS(cloud_counter_phase2,k+5) = Cloud(istep).fC;
            VLOSS(cloud_counter_phase2,k+6) = Cloud(istep).a;
            VLOSS(cloud_counter_phase2,k+7) = Cloud(istep).b;
            for k2 = 1:k
                VLOSS(cloud_counter_phase2,k2+4) = Cloud(idx_alt).Sed(k2) - Cloud(istep).Sed(k2);
            end
            fprintf(fid2, formatSpecVL, VLOSS(cloud_counter_phase2,1), VLOSS(cloud_counter_phase2,2), VLOSS(cloud_counter_phase2,3), VLOSS(cloud_counter_phase2,4), VLOSS(cloud_counter_phase2,5), VLOSS(cloud_counter_phase2,6), VLOSS(cloud_counter_phase2,7), VLOSS(cloud_counter_phase2,8), VLOSS(cloud_counter_phase2,9), VLOSS(cloud_counter_phase2,10)); % print current vloss results in output file
            cloud_counter_phase2 = cloud_counter_phase2 + 1;                    % increase the counter
            idx_alt = istep;
        end

        %% Updating Live GUI parameters for Dynamic Collapse
        if t_sim_istep == 1 && mod(dt*(istep-1),dt_sim) == 0 ...
                || t_sim_istep == 0 && mod(istep,25) == 0                       % update Live parameters in the GUI for Phase 1 for every xth (25) step

            % Set axis limits of live plot in GUI
            app.UIAxes.YLim = [-app.WaterdepthEditField.Value 0];               % set y limit to be between water surface and the water depth
            if max(max(x_half_sphere)) > app.WaterdepthEditField.Value          % extend the x axis by 50 meters every time the maximum displayable width is reached
                xlim_up = 50 * ceil(max(max(x_half_sphere))/50);                % round to the nearest multiple of 50
                app.UIAxes.XLim = [0 xlim_up];                                  % update the axis limits in the GUI
            else
                app.UIAxes.XLim = [0 app.WaterdepthEditField.Value];
            end
            % update GUI parameters
            app.TimeEditField.Value = round(Cloud(istep).t,1);                  % update current simulation time
            app.CloudwidthEditField.Value = round(Cloud(istep).b*2,1);          % update cloud width
            app.CloudheightEditField.Value = round(Cloud(istep).a,1);           % update cloud radius
            app.SettlingEditField.Value = sum(sum(VLOSS(cloud_counter_phase1:cloud_counter_phase2-1,5:(k+4)))); % update settling volume
            drawnow                                                             % show updated parameters in GUI
            % plot cloud contour
            [x_sphere,y_sphere,z_sphere] = sphere;                              % create a sphere
            r_ph2_plot = Cloud(istep).b;                                        % assign radius of the sphere (ellipsoid)
            r2_ph2_plot = Cloud(istep).a;                                       % assign radius of the sphere (ellipsoid)
            x_half_sphere = x_sphere(1:11,[1 11]).* r_ph2_plot;                 % get x-coordinates of the cloud sphere (ellipsoid)
            y_half_sphere = y_sphere(1:11,[1 11]).* r_ph2_plot;                 % get y-coordinates of the cloud sphere (ellipsoid)
            z_half_sphere = z_sphere(1:11,[1 11]).* r2_ph2_plot;                % get z-coordinates of the cloud sphere (ellipsoid)
            z_top_ph2_plot = Cloud(istep).a - 3/8 * r2_ph2_plot;                % get centre of mass coordinate
            z_half_ellipse = z_half_sphere*(-1) - app.WaterdepthEditField.Value; % get z-coordinates of the upper half ellipsoid
            hold(app.UIAxes,"off");                                             % axes setting
            plot(app.UIAxes,x_half_sphere,z_half_ellipse,'Color','k')           % plot cloud contour
            drawnow                                                             % show updated plot in GUI
        end


        %% Check criteria for termination of Phase 2
        % Calculate cloud velocity
        u_Cloud = sqrt(KE_P2(E_counter-1)*2/(Cloud(istep).m));                  % calculate velocity of the cloud by means of the kinetic energy and cloud mass

        % Calculate ambient velocity
        u_amb_avg = sqrt(sum(amb.^2));                                          % Average ambient velocity
        u_amb_crit = u_amb_crit;                                                % turbulent velocity that is assumed to be present in coastal waters even during slack tide
        u_amb = max(u_amb_avg,u_amb_crit);                                      % get maximum ambient velocity component

        if u_Cloud <= u_amb                                                     % exit criteria - compare velocities
            t_final = t(istep);                                                 % store final simulation time

            feedText = [feedText; 'End of Phase 2'];
            app.commandfeedTextArea.Value = feedText;                           % give feedback to GUI command feed
            disp('"End of Phase 2"')                                            % give feedback to command window
            t(end);                                                             % give feedback to command window (end time of Dynamic Collapse)

            fprintf(fid, ' -------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n'); % seperation line in output file
            fprintf(fid, ' End of Phase 2 after %5.2f seconds \n',t_final);     % write in output file
            fprintf(fid, ' -------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n'); % seperation line in output file
            fclose(fid);

            fprintf(fid2, ' -------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n'); % seperation line in output file
            fprintf(fid2, ' End of Phase 2 after %5.2f seconds \n',t_final);    % write in output file
            fprintf(fid2, ' -------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n'); % seperation line in output file
            fclose(fid2);

            break
        end

        %% Calculation of the gradients with the numerical solver
        [Cloud_Input] = prover_m_rk4(Cloud_Input,b_impact,Hopper,Sediment,Para,d,dt,qe,amb,istep,z_i,iz,rho_amb); % enter numerical solver (RK4)

        %% prepare next time step
        t(istep+1)=istep*dt;                                                    % time for next time step
        istep=istep+1;                                                          % increase the time step counter

    end

    %% Save variables
    save(char(strcat(file_info{1},'/output/',name,'_dump',num2str(i_dump),'.mat')),'Cloud','VLOSS');      % save the cloud and stripped/settled material

end

end