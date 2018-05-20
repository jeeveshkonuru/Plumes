function [ S ] = thermal_adjustment( S )
if S.thermal == 0
    S.rho_wt = S.rho_w;
    return
end

if S.constprof == 0
    S.Re=(2.*S.flowrate.*S.rho0)./(S.mu.*S.bi.*pi); %Reynolds number in pipe
    f=0.184*(S.Re).^(-1/5); %friction factor, smooth pipe
    d=2.*S.bi; %Pipe Diameter
    S.Pr = S.cp*S.mu/S.kw; %Prandtl number
    S.Nu= f/8*(S.Re-1000)*S.Pr/(1+12.7*sqrt(f/8)*(S.Pr^(2/3)-1)); %Nusselt number
    h = S.Nu*S.kw/d;
    mdot = S.flowrate*S.rho_w;
    S.KK = pi*d/(mdot*S.cp)*h; %Total coefficient for the temperature relation
    
    dT = @(z,T) S.KK.*(S.T_a(z) - T);
    updT = @(z,T) -dT(z,T);
    
    if S.surfacetemp == 1
        Tsurf = S.T_a(S.zmin);
    else
        [zout_up, Tout_up] = ode45(updT,[S.zmax,S.zmin+1],S.T_a(S.zmax));
        Tsurf= Tout_up(end);
    end
    [zout, Tout] = ode45(dT,[S.zmin+1 , S.zr],Tsurf);
    [zout_1, Tout_1] = ode45(dT,[S.zmin+1 , S.zr],S.T_a(S.zmin));
    
    S.Tr = Tout(end);
    S.Tru = Tout_1(end);
    
    %Below just assumes comes from top or bottom; Doesn't take into accoutn
    %given density. NEED TO FIX
    
    SA_ud=gsw_SA_from_SP(S.psal_a(S.zmin),S.zr,S.long,S.lat);
    CT_ud = gsw_CT_from_t(SA_ud,S.Tru,S.zr);
    S.rho_wtu = gsw_rho(SA_ud,CT_ud,0);
    
    SA_u=gsw_SA_from_SP(S.psal_a(S.zmax),S.zr,S.long,S.lat);
    CT_u = gsw_CT_from_t(SA_u,S.Tr,S.zr);
    S.rho_wt = gsw_rho(SA_u,CT_u,0);
    
    
elseif S.constprof == 1
    error = 'thermal adjustment not implemented yet not implemented yet'
else
    error = 1


end

