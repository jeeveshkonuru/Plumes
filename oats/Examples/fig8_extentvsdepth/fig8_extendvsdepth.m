%Andrew Rzeznik
%Code for running plume extent vs volume flux, all other factors held
%constant

% clear
% clc
% close all

s1=param_file_to_struct('fig8_extentvsdepth_s');
s2=param_file_to_struct('fig8_extentvsdepth_clim_s');
s3=param_file_to_struct('non_thermal_top_s');
S=plumes_main(s1);
Sclim=plumes_main(s2);
St=plumes_main(s3);
toldif=10;

Ls0=length(St);
Ls=length(S);
depth=[];
depthclim=[];
zrel=[];
zrelclim=[];

zr=zeros(Ls0,1);
dtop=zr;

for i=1:Ls0
    zr(i)=St{i}.zr;
    dtop(i) = St{i}.Z;
end
zrt=zr;

i=1;
while i < Ls0
    if real(dtop(i)) > 400
        dtop=[dtop(1:i-1); dtop(i+1:end)];
        zrt=[zrt(1:i-1); zrt(i+1:end)];
        Ls0 = Ls0-1;
    else
        i=i+1;
    end
end

for i = 1:Ls
    if S{i}.Z < 300
        zrel = [zrel S{i}.zr];
        depth = [depth S{i}.Z];
    end
    if Sclim{i}.Z < 350
        zrelclim = [zrelclim Sclim{i}.zr];
        depthclim = [depthclim Sclim{i}.Z];
    end
end




