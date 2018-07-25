function idSocial_auxiliaries_displayPlotBox(ax)

if nargin <1 || isempty(ax)
    ax=gca;
end

% if verLessThan('matlab','8.5.0')
%     box(ax,'off')
%     %     set(ax,'Color','none')
%     fakeax=axes('Position',get(ax,'Position') ,'box','on','xtick',[],'ytick',[],'LineWidth',get(ax,'LineWidth')+10,'Color','none');
%     uistack(fakeax,'top')
% else
    %     ax.Box = 'off';
%     fakeax=axes('Units',get(ax,'Units'),'Position',get(ax,'Position') + [1 1 0 0],'box','on','xtick',[],'ytick',[],'LineWidth',get(ax,'LineWidth'),'Color','none');
%     box(fakeax,'on')
%     uistack(fakeax,'top')
    box(ax,'off')
%     set(ax,'Color','none')
    
    xlim = get(ax,'XLim');
    ylim = get(ax,'YLim');
    axes(ax);
    hold on
    plot(xlim,[ylim(2) ylim(2)],'-','Color','k','LineWidth',get(ax,'LineWidth'));
    plot([xlim(2) xlim(2)],ylim,'-','Color','k','LineWidth',get(ax,'LineWidth'));

    set(ax,'YLim',ylim)
    set(ax,'XLim',xlim)
% end


