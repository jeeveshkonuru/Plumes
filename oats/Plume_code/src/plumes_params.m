function S=plumes_params(sout)
%Andrew Rzenzik 8/9/2016
%Function which constructs the runs for the plume based on a given
%parameter file


%Brings in the initial struct with all variables
%eval(filename) NO LONGER USED

%% Construct Cell array of runs
%Unwinds the struct input data to create a cell array with one run struct
%in each entry

%Use the runmaker recursive function to unwind the struct into a list of
%structs, each of which represent the input data for a single run.
S=run_maker({sout});


%% Read the data out of the given file for stratification
%This could be heavily optimized to save time outside the loop by checking
%if the filename changes 
numruns=length(S);

%The loop below calculated the ambient potential density and other initial
%conditions for the parameters for each run

%% 
%Below code is **TEMPORARY, builds a few currently non-iteratable parameters.
%Ideally this would later take filenames in the struct, iterating over each
%struct in the cell array, and replace them with the associated data
% 
% numruns=length(S);

% casenumber=2;%**TEMPORARY Casenumber sets the background stratification
% [zp,rhop]=rhotestc(casenumber); %**TEMPORARY SOLUTION will need better code for this
% 
%% Non-dimensionalized
%Below code goes through and create the non-dimensionalized parameters, or
%other stuff

for i=1:numruns
    if S{i}.climatology == 1
        %%Load the data from the argo climate file
        load(S{i}.profile);
        pres_raw=pressure_mat;
        rho=rho_mat;
        rhotheta=rhotheta_mat;
        temp_raw=temp_mat;
        psal_raw=SA_mat;
        longitude = long_mat;
        latitude = lat_mat;
        %%Finish density construction
        
        drhothetadz = diff(rhotheta)./diff(pres_raw); %Gradient calculated at midpoints
        pres_mid = (pres_raw(2:end) + pres_raw(1:(end-1)))./2; %Pressure at the midpoints of the grid
        %% Constructing temperature profiles and height over which they are valid,
        %and initial and end parameters.
        S{i}.zmin = pres_raw(1);
        S{i}.zmax = pres_raw(end);
        S{i}.T_a =@(z) interp1q(pres_raw,temp_raw,z);
        S{i}.psal_a = @(z) interp1q(pres_raw,psal_raw,z);
        S{i}.long = longitude;
        S{i}.lat = latitude;
    
        S{i}.rhoa=@(z) interp1q(pres_raw,rho,z); %Linear piecewise density function via interpolation
        S{i}.rhoatheta = @(z) interp1q(pres_raw,rhotheta,z); %Same for potential density
        S{i}.N2 = @(z) S{i}.g./S{i}.rho0.*interp1q(pres_mid,drhothetadz,z); %Linear piecewise bouyancy frequency
        S{i}.c = 1484; %Speed of sounds; currently is constant, but is only used for initial density
        S{i}  = thermal_adjustment( S{i} );
        S{i}.rhoi=(S{i}.m_s+S{i}.m_f+S{i}.rho_w.*S{i}.flowrate)./(S{i}.m_s./S{i}.rho_s+S{i}.m_f./S{i}.rho_f+S{i}.flowrate.*S{i}.rho_w./S{i}.rho_wt); %Calculating the initial density at the release depth
        S{i}.ui = S{i}.flowrate./(S{i}.bi.^2.*pi); % Calculates initial velocity parameter based on flow rate, and radius. ** May change based on the profile being modelled.
        S{i}.Qi = S{i}.flowrate;
        S{i}.Mi = S{i}.Qi.*S{i}.ui;
        S{i}.Fi = S{i}.Qi.*S{i}.g.*(S{i}.rhoi-S{i}.rhoatheta(S{i}.zr))./S{i}.rho0;
    elseif S{i}.constprof == 1
        S{i}.rhoatheta = @(z) S{i}.rho_water+S{i}.rho0./S{i}.g.*S{i}.N2(z).*z; %Same for potential density
        S{i}.rhoa = @(z) S{i}.rhoatheta(z) + S{i}.rho0.*S{i}.g.*z./S{i}.c.^2;
        S{i}.rhoi=S{i}.rhoatheta(S{i}.zr)+S{i}.rho_waste; %Calculating the initial density at the release depth
%        S{i}.zmax=max(pres_raw); %Sets the maximum depth of the plume run in the integrator
        S{i}.ui = S{i}.flowrate./(S{i}.bi.^2.*pi); % Calculates initial velocity parameter based on flow rate, and radius. ** May change based on the profile being modelled.
        S{i}.Qi = S{i}.flowrate;
        S{i}.Mi = S{i}.Qi.*S{i}.ui;
        S{i}.Fi = S{i}.Qi.*S{i}.g.*(S{i}.rhoi-S{i}.rhoa(S{i}.zr))./S{i}.rho0;
    elseif S{i}.constprof == 0 
    
        %Code which reads and argo profile
        [latitude,longitude,positioning_system,position_qc,juld_location,...
        juld,juldqc,pressure,pres_raw,pres_qc_raw,pres_qc,pressure_error,...
        pres_units, pres_fillvalue,temperature,temp_raw,temp_qc_raw,temp_qc,...
        temperature_error,temp_units,temp_fillvalue...
        salinity,psal_raw,psal_qc_raw,psal_qc,salinity_error,psal_units,psal_fillvalue...
        platform_number,cycle_number,data_type,format_version,project_name,...
        pi_name,direction,data_centre,data_state_indicator,data_mode,...
        dc_reference,platform_type,wmo_inst_type,station_parameters,...
        reference_date_time,vertical_sampling,config_mission_number]=...
        argo_profile_read_matlab2008bplus(S{i}.profile);
    
        %below code calculates the potential density from a given ARGO profile,
        %along with a gradient which is continuous and piecewise linear.
        SA=gsw_SA_from_SP(psal_raw,pres_raw,longitude,latitude);
        CT = gsw_CT_from_t(SA,temp_raw,pres_raw);
        rho = gsw_rho(SA,CT,pres_raw); %Calculate density from Salinity, temp, pressure
        rhotheta = gsw_rho(SA,CT,0); % Potential density by setting pressure to zero
        drhothetadz = diff(rhotheta)./diff(pres_raw); %Gradient calculated at midpoints
        pres_mid = (pres_raw(2:end) + pres_raw(1:(end-1)))./2; %Pressure at the midpoints of the grid
        %% Constructing temperature profiles and height over which they are valid,
        %and initial and end parameters.
        S{i}.zmin = pres_raw(1);
        S{i}.zmax = pres_raw(end);
        S{i}.T_a =@(z) interp1q(pres_raw,temp_raw,z);
        S{i}.psal_a = @(z) interp1q(pres_raw,psal_raw,z);
        S{i}.long = longitude;
        S{i}.lat = latitude;
    
        S{i}.rhoa=@(z) interp1q(pres_raw,rho,z); %Linear piecewise density function via interpolation
        S{i}.rhoatheta = @(z) interp1q(pres_raw,rhotheta,z); %Same for potential density
        S{i}.N2 = @(z) S{i}.g./S{i}.rho0.*interp1q(pres_mid,drhothetadz,z); %Linear piecewise bouyancy frequency
        S{i}.c = 1484; %Speed of sounds; currently is constant, but is only used for initial density
        S{i}  = thermal_adjustment( S{i} );
        S{i}.rhoi=(S{i}.m_s+S{i}.m_f+S{i}.rho_w.*S{i}.flowrate)./(S{i}.m_s./S{i}.rho_s+S{i}.m_f./S{i}.rho_f+S{i}.flowrate.*S{i}.rho_w./S{i}.rho_wt); %Calculating the initial density at the release depth
        S{i}.ui = S{i}.flowrate./(S{i}.bi.^2.*pi); % Calculates initial velocity parameter based on flow rate, and radius. ** May change based on the profile being modelled.
        S{i}.Qi = S{i}.flowrate;
        S{i}.Mi = S{i}.Qi.*S{i}.ui;
        S{i}.Fi = S{i}.Qi.*S{i}.g.*(S{i}.rhoi-S{i}.rhoatheta(S{i}.zr))./S{i}.rho0;
        
    else
        error = 1
    end
end

% for i=1:numruns
   
%     elseif S{i}.standardprof == 1
%         S{i}.c = 1484; %Speed of sound in water, m/s **This currently is just approximated as a constant, but should be able to be calculated based on position and SEAWATER TABLES
%         S{i}.rhoi=S{i}.rho_water+S{i}.rho_waste+S{i}.rho0.*S{i}.g.*S{i}.zr./S{i}.c.^2; %Initial density of plume release **USES A CONSTANT speed c
%         S{i}.rhoa=@(z) interp1q(zp,rhop,z); %The ambient density curve
%         S{i}.zmax=max(zp); %Sets the maximum depth of the plume run in the integrator
%         S{i}.G=S{i}.g.^2./S{i}.c.^2;
%         S{i}.F0=1/S{i}.rho0;
%         S{i}.w0=S{i}.flowrate/pi*2.^(-7/8)*pi^(3/4)*S{i}.entr^(-1/2)*S{i}.F0^(-3/4)*S{i}.G^(5/8); %Initial non-dimensionalized horizontal paramter
%         S{i}.v0=S{i}.flowrate/(pi*S{i}.bi)*2^(-3/4)*pi^(1/2)*S{i}.F0^(-1/2)*S{i}.G^(1/4); %Initial non-dimentionalized velocity
%         S{i}.f=-1; %Initial non-dimensionalized Bouyancy
%         S{i}.N20 = S{i}.g.*0.0051./S{i}.rho0;
%         S{i}.G = S{i}.N20 - (S{i}.g./S{i}.c).^2; %Plus or minus?
%     end
% end

end