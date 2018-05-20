%parameters for figure and panel size
%Thanks Patrick Martineau

%% Run the code to generate/load the data here
%Setting up the initial data
clear
clc
close all
fig7_standardplume




%% Calling parameter functions and setting up initial plot
plot_params_fig7;
plot_setup;
 
%% loop to create axes and plots
for i=1:subplotsx
    for ii=1:subplotsy
        ax=axes('position',sub_pos{i,ii},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
        set(ax, 'Ydir', 'reverse') %OPTIONAL
        hold on
        
        %Code for each subplot ----------------------------------
        if i == 1
            %Radius
            plot(S{1}.Q./sqrt((S{1}.M)),S{1}.Z-1000,'k')
            axis([0 60 0 220])
            grid
            xlabel('b (m)')
            ylabel('z (m)')

        elseif i == 2
            %mean turbulentVelocity
            plot(S{1}.M./(S{1}.Q),S{1}.Z-1000,'k')
            axis([0 1 0 220])
            grid
            xlabel('u (m/s)')
            set(ax,'yticklabel',[])
            
        elseif i == 3
            %Potential Density anomily
            plot(S{1}.rho0./S{1}.g.*S{1}.F./(S{1}.Q),S{1}.Z-1000,'k')
            axis([-0.5 1 0 220])
            xlabel('\Delta\rho_\theta (kg/m^3)')
            set(ax,'yticklabel',[])
            grid

        end
        %-------------------------------------------------------------------
        %set(ax,'yticklabel',[])
        %Labels the figures
        text(-.07*(ax.XLim(2)-ax.XLim(1))+ax.XLim(1),-0.07*(ax.YLim(2)-ax.YLim(1))+ax.YLim(1),['(',char(i+96),')'],...
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