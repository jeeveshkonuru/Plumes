%Andrew Rzeznik
%Code for running plume extent vs volume flux, all other factors held
%constant

clear
clc
close all

param_file = 'fig7_standardplume_s';
sout = param_file_to_struct(param_file);
S=plumes_main(sout);

