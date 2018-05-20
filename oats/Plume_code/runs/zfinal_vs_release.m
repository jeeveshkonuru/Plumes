%Andrew Rzeznik
%Code for running plume extent vs volume flux, all other factors held
%constant

clear
clc
close all

param_file = 'zfinal_vs_release_s';
S=plumes_main(param_file);

toldif=10;

Ls0=length(S);
a0=zeros(Ls0,1);
b0=zeros(Ls0,1);
for i=1:Ls0
    a0(i)=S{i}.Z;
    b0(i)=S{i}.zr;
end



%Checking for erros with jumps in dif eq
i=1;
%Below code doesn't currently handle edge cases FIX THIS
while i < Ls0
    if a0(i)-a0(i+1) > toldif
        a0=[a0(1:i-1); a0(i+1:end)];
        b0=[b0(1:i-1); b0(i+1:end)];
        Ls0 = Ls0-1;
    else
        i=i+1;
    end
end

%post processing for this case the correct values
% pikind = real(a0) > 0.08;
% a0=a0(pikind);
% b0=b0(pikind);


figure
hold1=axes;
plot(a0,b0);
title('Final plume vertical extent based on release depth')
xlabel('Plume vertical extent (m)')
ylabel('Release depth (m)')
grid
set(hold1, 'Ydir', 'reverse')
