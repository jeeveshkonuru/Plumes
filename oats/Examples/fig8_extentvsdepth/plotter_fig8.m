%parameters for figure and panel size
%Thanks Patrick Martineau

%% Run the code to generate/load the data here
%Setting up the initial data
clear
clc
close all
fig8_extendvsdepth




%% Calling parameter functions and setting up initial plot
plot_params_fig8;
plot_setup;
 
%% loop to create axes and plots
for i=1:subplotsx
    for ii=1:subplotsy
        ax=axes('position',sub_pos{i,ii},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
        %set(ax, 'Ydir', 'reverse') %OPTIONAL
        hold on
        
        %Code for each subplot ----------------------------------
        if i == 1
            % extent vs depth, 
            plot(zrel,depth,'k')
            plot(zrelclim,depthclim,'k--')
            plot(zrt,dtop,'k:')
            xlabel('z_i (m)')
            ylabel('L(m)')
            legend({'Single ARGO stratification','Average summer stratification', 'Single stratification surface temperature'},'Location','southeast')
            grid
        end
        %-------------------------------------------------------------------
        %set(ax,'yticklabel',[])
        %Labels the figures
        %text(-.07*(ax.XLim(2)-ax.XLim(1))+ax.XLim(1),1.07*(ax.YLim(2)-ax.YLim(1))+ax.YLim(1),['(',char(i+96),')'],...
        %    'color','k',...
        %    'fontw','b')
        
        
        
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