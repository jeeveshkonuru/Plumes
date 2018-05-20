        % FIGURE POSTPROCESSING
        %Labels the figures

        
        
        
        % text(1,1,['(',char(ii+96),') '],...
        %     'HorizontalAlignment','right',...
        %     'VerticalAlignment','top',...
        %     'color','b',...
        %     'fontw','b')
        
        
        
        %Sets font for current plot
        set(gca,'FontSize',fontsize)
        set(findall(gcf,'type','text'),'FontSize',fontsize)