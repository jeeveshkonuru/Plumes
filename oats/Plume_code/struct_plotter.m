Ls=length(S);
Tvec=zeros(Ls,1);
Tavec=Tvec;
zvec=zeros(Ls,1);
for i=1:Ls
    Tvec(i)=S{i}.rho_wt;
    Tavec(i)=S{i}.rhoatheta(S{i}.zr);
    zvec(i)= S{i}.zr;
end
figure
hold1=axes;
plot(Tvec-Tavec,zvec);
set(hold1, 'Ydir', 'reverse')


    