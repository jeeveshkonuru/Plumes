%parameters for figure and panel size
%Thanks Patrick Martineau

%% Run the code to generate/load the data here
%Setting up the initial data
clear
clc
close all





%% Calling parameter functions and setting up initial plot
plot_params_fig2;
plot_setup;
 
%% loop to create axes and plots
for i=1:subplotsx
    for ii=1:subplotsy
        ax=axes('position',sub_pos{i,ii},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
        set(ax, 'Ydir', 'reverse') %OPTIONAL
        hold on
        
        %Code for each subplot ----------------------------------
        if i == 1
            % Ambient Potential Density
            plot(rhop,zp)
            plot(rho,pres_raw,'--')
            xlabel('\rho_{\theta} (kg/m^3)')
            ylabel('z (m)')
            grid
        elseif i == 2
            
            %Ambient Potential Density Gradient
            plot(drhoadz(zp(2:end-1),zp,rhop),zp(2:end-1))
            plot(drhoadz(pres_raw(2:end-1),pres_raw,rho),pres_raw(2:end-1),'--')
            
            grid
            xlabel('\partial\rho_\theta/\partial z')
            
            set(ax,'yticklabel',[])
            
            text(-.05*(ax.XLim(2)-ax.XLim(1))+ax.XLim(1),-0.05*ax.YLim(2),['(',char(i+96),')'],...
                'color','k',...
                'fontw','b')
            %The inlaid plot
            axis([0 0.06 0 2000])
            h3 = axes('Position',[ax.Position(1)+ax.Position(3)*0.27 ax.Position(2)+ax.Position(4).*0.12  ax.Position(3).*2./3 ax.Position(4)*0.625],'Fontsize',fontsize);
            box on
            plot(drhoadz(zp(40:end-1),zp,rhop),zp(40:end-1))
            hold on
            plot(drhoadz(pres_raw(50:end-1),pres_raw,rho),pres_raw(50:end-1),'--')
            set(h3, 'Ydir', 'reverse')
            grid
            %xlabel('\partial\rho_\theta/\partial z')
            %ylabel('z')
            
            param_file = 'Theta_calc_clim_s';
            S=plumes_main(param_file);
            Sspec=plumes_main('Theta_calc_s');
            
            toldif=10;
            
            Ls0=length(S);
            Ls0av=length(Sspec);
            zwrite=(10:10:1900)';
            
            a0=S{1}.N2(zwrite).*S{1}.c.^2./S{1}.g.^2+1;%theta
            b0=zwrite;
            
            a0spec=Sspec{1}.N2(zwrite).*Sspec{1}.c.^2./Sspec{1}.g.^2+1;%theta
            
            
        elseif i == 3
            %Compressibility number
            plot(a0,b0);
            plot(a0spec,b0,'--')
            xlabel('\Theta')
            grid
            legend('Average','Single-profile','Location','southeast')
            
        end
        %-------------------------------------------------------------------
        %set(ax,'yticklabel',[])
        %Labels the figures
        text(-.07*(ax.XLim(2)-ax.XLim(1))+ax.XLim(1),1.07*(ax.YLim(2)-ax.YLim(1))+ax.YLim(1),['(',char(i+96),')'],...
            'color','k',...
            'fontw','b')
        
        
        
        % text(1,1,['(',char(ii+96),') '],...
        %     'HorizontalAlignment','right',...
        %     'VerticalAlignment','top',...
        %     'color','b',...
        %     'fontw','b')
        
        
        
        %Sets font for current plot
        plot_postprocess
    end
end

%% Saving eps with matlab and then producing pdf and png with system commands
plot_output;