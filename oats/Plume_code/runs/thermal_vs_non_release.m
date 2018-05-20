%Andrew Rzeznik
%Code for running plume extent vs volume flux, all other factors held
%constant

clear
clc
close all

S=plumes_main(param_file);
St=plumes_main(param_file);
Sb=plumes_main(param_file);

toldif=10;

Ls0=length(S);
zwrite=(10:10:1900)';

a0=S{1}.N2(zwrite).*S{1}.c.^2./S{1}.g.^2+1;%theta
b0=zwrite;%depth
% for i=1:Ls0
%     a0(i)=S{i}.Z;
%     b0(i)=S{i}.zr;
% end



%Checking for erros with jumps in dif eq
i=1;
%Below code doesn't currently handle edge cases FIX THIS
% while i < Ls0
%     if a0(i)-a0(i+1) > toldif
%         a0=[a0(1:i-1); a0(i+1:end)];
%         b0=[b0(1:i-1); b0(i+1:end)];
%         Ls0 = Ls0-1;
%     else
%         i=i+1;
%     end
% end

%post processing for this case the correct values
% pikind = real(a0) > 0.08;
% a0=a0(pikind);
% b0=b0(pikind);


figure
hold1=axes;
plot(a0,b0);
title('Compressibility number for a given density profile')
xlabel('\Theta')
ylabel('Depth (m)')
grid
set(hold1, 'Ydir', 'reverse')
