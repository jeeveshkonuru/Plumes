function [z,Q,M,F,x] = plumes_run(sout)
%PLUME_RUN 8/9/2016
%Runs the dimensional mass-matrix form of the plumes function given the
%initial data, outputing dimensional parameters for the solution
%CHANGES: UA actually a fcn of v4
VEn = @(v) ((sout.entr.*abs(v(2))).^sout.entr_num+(sout.entrb.*sout.Ua(v(4),sout.U_para,sout.U_parb).*v(1)).^sout.entr_num).^(1./sout.entr_num);
%2.*abs(v(2).*v(1).*En(v)./(v(2).^2+v(1).^2*sout.Ua.^2).^(3/4))
%2.*sout.entr.*sqrt(v(2))
f=@(z,v) [2.*abs(v(2)).*VEn(v)./(v(1).^2.*sout.Ua(v(4),sout.U_para,sout.U_parb).^2+v(2).^2).^(3/4);v(1).*v(3)./sqrt(v(2).^2+v(1).^2*sout.Ua(v(4),sout.U_para,sout.U_parb).^2);-sout.N2(v(4)).*v(1).*v(2)./sqrt(v(2).^2+v(1).^2*sout.Ua(v(4),sout.U_para,sout.U_parb).^2);v(2)./sqrt(v(2).^2+v(1).^2*sout.Ua(v(4),sout.U_para,sout.U_parb).^2);v(1).^2*sout.Ua(v(4),sout.U_para,sout.U_parb).^2./sqrt(v(2).^2+v(1).^2*sout.Ua(v(4),sout.U_para,sout.U_parb).^2)]; %NEED to add two extra equations for position-added at end

%options = odeset('Mass',M); % Sets the basic options setup for the ODE

% runfun=@(z,v,flag) mass_ode(z,v,flag,sout);
[zout,vout] = ode45(f,[0 5000],[sout.Qi,sout.Mi,sout.Fi,sout.zr,0]'); %Changed with hardcoded integration parameter

Qt=vout(:,1);
Mt=vout(:,2);
Ft=vout(:,3);
Zt=vout(:,4);
Xt=vout(:,5);

Zt=Zt(imag(Qt) == 0);
Qt=Qt(imag(Qt) == 0);
Mt=Mt(imag(Qt) == 0);
Ft=Ft(imag(Qt) == 0);
Xt=Xt(imag(Qt) == 0);

%Output code
if sout.fullprof==1
    z=Zt;
    Q=Qt;
    M=Mt;
    F=Ft;
    x=Xt;
else %Not Done Yet, FIX (interpolation issue)
    [zabsmax,ind] = max(Zt);
    z=zabsmax-sout.zr;
    Q=Qt(ind);
    M=Mt(ind);
    F=Ft(ind);
    x=Xt(ind);
end
end


