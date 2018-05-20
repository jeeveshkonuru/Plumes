function S=plumes_main(Sin)
%Andrew Rzeznik
%Mining Waste Plumes Code
%8/9/2016
%A Test script that should test all functionality of the plumes software.
%Also used as the main jumping off point in implementation of the software
%
%To run this code, 

%% Clearing Code and setting parameters

%THIS SECTION SHOULD BE COMMENTED OUT IF THE SCRIPT IS BEING USED IN A
%LARGER SCRIPT TO GENERATE DATA

% clear
% clc
% param_file = 'testvals';

%% Main - Runs the simulations

%Below function constructs the run parameters from the param_file

S=plumes_params(Sin);

%Below code makes all the runs and outputs the results based on the desired
%output of the param file
[Z,Q,M,F,X]=cellfun(@plumes_makeprof,S,'UniformOutput',0);
Zd=zeros(1,length(Z));
Bd=Zd;
Ud=Zd;
RHOd=Zd;
znd=Zd;
vnd=Zd;
wnd=Zd;
zndi=Zd;
vndi=Zd;
wndi=Zd;


for i = 1:length(S)
    if S{i}.nondim == 1 && S{i}.fullprof==0 %NONDIM NOT IMPLEMENTED
        S{i}.Zd = 2.^(-7/8).*pi.^(-1/4).*S{i}.entr.^(-1/2).*S{i}.F0.^(1/4).*S{i}.G.^(-3/8).*Z(i);
        S{i}.V = 2.^(3/4).*pi.^(-1/2).*S{i}.F0.^(1/2).*S{i}.G.^(-1/4).*U(i);
        S{i}.W = 2.^(7/8).*pi.^(-3/4).*S{i}.entr.^(1/2).*S{i}.F0.^(3/4).*S{i}.G.^(-5/8).*B(i);
        S{i}.F = 2.*pi.^(-1).*S{i}.F0.*RHO(i);
        Zd(i) = S{i}.Zd;
        V(i) = S{i}.V;
        W(i) = S{i}.W;
        F(i) = S{i}.F;
        S{i}.Bd = S{i}.W./S{i}.V;
        S{i}.Ud = S{i}.V.^2./S{i}.W;
        S{i}.RHOd = 2./pi.*S{i}.F.*(S{i}.rho0./(S{i}.g.*S{i}.W))+S{i}.rhoa(S{i}.Zd.*S{i}.finalfrac);
        Bd(i) = S{i}.Bd;
        Ud(i) = S{i}.Ud;
        RHOd(i) = S{i}.RHOd;
    elseif S{i}.nondim == 0
        S{i}.Z=Z{i};
        S{i}.Q=Q{i};
        S{i}.M=M{i};
        S{i}.F = F{i};
        S{i}.X = X{i};
%         S{i}.wnd = S{i}.W./(2.^(7/8).*pi.^(-3/4).*S{i}.entr.^(1/2).*S{i}.F0.^(3/4).*S{i}.G.^(-5/8));
%         S{i}.vnd = S{i}.V./(2.^(3/4).*pi.^(-1/2).*S{i}.F0.^(1/2).*S{i}.G.^(-1/4));
%         S{i}.znd = Z(i)./(2.^(-7/8).*pi.^(-1/4).*S{i}.entr.^(-1/2).*S{i}.F0.^(1/4).*S{i}.G.^(-3/8));
%         znd(i) = S{i}.znd ;
%         wnd(i) = S{i}.wnd;
%         vnd(i) = S{i}.vnd;
        %Initial parameters
        S{i}.Mi=S{i}.bi.^2.*S{i}.ui.^2;
%       S{i}.Qi=S{i}.bi.^2.*S{i}.ui;
%         S{i}.Fi = S{i}.Qi.*S{i}.g.*(S{i}.rhoa(Z(i))-S{i}.rhoi)./S{i}.rho0;
%         S{i}.wndi = S{i}.Wi./(2.^(7/8).*pi.^(-3/4).*S{i}.entr.^(1/2).*S{i}.F0.^(3/4).*S{i}.G.^(-5/8));
%         S{i}.vndi = S{i}.Vi./(2.^(3/4).*pi.^(-1/2).*S{i}.F0.^(1/2).*S{i}.G.^(-1/4));
%         S{i}.zndi = S{i}.zr./(2.^(-7/8).*pi.^(-1/4).*S{i}.entr.^(-1/2).*S{i}.F0.^(1/4).*S{i}.G.^(-3/8));
%         zndi(i) = S{i}.zndi ;
%         wndi(i) = S{i}.wndi;
%         vndi(i) = S{i}.vndi;
    end
end
%% Plotting Code
% fig=figure;
% plot(zr,zfinal)
% title('Plume vertical extent vs. release depth')
% xlabel('Release depth')
% ylabel('Plume vertical extent (m)')
% grid
% 
% filename='z_vs_zr';
% savefig(fig,filename)
% saveas(fig,filename,'eps')
% saveas(fig,filename,'png')
% 
% fig=figure;
% plot(zr,2.*wfinal)
% title('Plume width vs. release depth')
% xlabel('Release depth')
% ylabel('Plume width (m)')
% grid
% 
% filename='w_vs_zr';
% savefig(fig,filename)
% saveas(fig,filename,'eps')
% saveas(fig,filename,'png')
% 
% save('release')