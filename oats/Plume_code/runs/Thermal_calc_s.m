%Fill out each of the following with cell arrays of the appropriate value
%type. MUST HAVE sout as variable

%ALL CHARACTER STRINGS MUST BE IN A CELL ARRAY IF MULTIPLE VALUES
%Currently the output loops over the lowest vector first, highest vector is
%highest loop. This should be codified

sout.rho0=1024; %Reference Density in kg/m^3
sout.g=9.81; % Local gravity acceleration m/s^2
sout.entr= 0.1; % Entrainment coefficient
sout.flowrate = 0.56; % Flowrate out of pipe, m^3/s
sout.zr = 10:10:1500; %Release depth of plume
sout.bi = 0.25 ;%Initial radial parameter; Corresponds to pipe radius in simple square profile case.
sout.rho_waste=10;%8.29; %Waste density, Kilogram/meter^3 THIS IS NOT CORRECT/USEFUL;REMOVE
sout.rho_f = 3334;%Density of manganese nodule fines, kg/m^3
sout.m_f = 1.54;%Mass flow rate of manganese nodule fines, in kg/s
sout.rho_s = 2750;%Density of sediment, in kg/m^3
sout.m_s = 3.1;% mass flow rate of sediment, in kg/s
sout.rho_w=1023.8;%Initial water density, at surface.;
sout.mu = 8.9e-4;%Viscosity of water
sout.kw = 600e-3;%Thermal conductivity of water
sout.cp= 4185.5; %Specific heat of water
sout.finalfrac = 0.9; % Because the plumes spread out to ininity and are not perfect near the end, this parameter, when just pulling the final depth or width, will decide what fraction down the plume to pull it from.
sout.fullprof=0;%Determines if the output should be the full profile or just a single value
sout.profile='R4901443_161.nc';%File for background profile. Currently assumes it is an ARGO single profile 
sout.constprof = 0; %use a constant profile instead of an ARGO profile
sout.N2 = @(z) 5e-6; %Bouyancy frequency squared for the constant stratification case
sout.c = 1484; %Speed of sounds; currently is constant, but is only used for initial density
sout.zmax = 1900; %Sets max depth of integrator
sout.nondim = 0; %runs non-dimensionalized code. NOT IMPLEMENTED CURRENTLY 
sout.thermal = 1; %Tells the code if it should use thermal calculations
sout.climatology = 0;%Tells the code if it should be using a climatology
sout.surfacetemp = 1;