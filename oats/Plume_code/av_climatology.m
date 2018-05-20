%Andrew Rzeznik
%12/2/1989
%File to load the variables for the climatology, to make a simpler .mat
%file

[pressure_mat,rho_mat,rhotheta_mat,long_mat,lat_mat,temp_mat,SA_mat] = argo_climatology_load();
save('argoclim_jul.mat');