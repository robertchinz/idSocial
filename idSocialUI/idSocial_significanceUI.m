function idSocial_significanceUI(lineH)

plotAx = get(lineH,'Parent');
plotfh = get(plotAx{1},'Parent');

mh = uimenu(plotfh,'Label','Statistics');
eh1 = uimenu(mh,'Label','Error bars');
eh11 = uimenu(eh1,'Label','Standard deviation','Callback',@addSTD);
eh12 = uimenu(eh1,'Label','SEM','Callback',@addSEM);
eh2 = uimenu(mh,'Label','Statistical Signifcance','Callback',@addStatSign);

linewidth = NaN(numel(lineH),1);
linecolor = NaN(numel(lineH),3);
xdataCell = cell(numel(lineH),1);
ydataCell = cell(numel(lineH),1);

for k = 1:numel(lineH)
    linewidth(k) = get(lineH(k),'LineWidth');
    linecolor(k,:) = get(lineH(k),'Color');
    set(lineH(k),'DeleteFcn',{@LineDelete,k});
    %     set(lineH(k), 'ButtonDownFcn', {@LineSelected,k})
    set(lineH(k), 'ButtonDownFcn', {@PointSelected,k})
    xdataCell{k} = get(lineH(k),'XData');
    ydataCell{k} = get(lineH(k),'YData');
end


marked_lines = zeros(1,numel(lineH));

marked_points = zeros(max(cellfun(@(x) numel(x),xdataCell)),numel(lineH));
markh = zeros(max(cellfun(@(x) numel(x),xdataCell)),numel(lineH));

no_possiblemarks = inf;


    function addSTD(src,event)
        errorbar(xticks,dat,nanstd(gui.plotResultsIndivData),'Color',get(ph,'Color'))
    end
    function addSEM(src,event)
        errorbar(xticks,dat,nanstd(gui.plotResultsIndivData)./sqrt(sum(~isnan(gui.plotResultsIndivData))),'Color',get(ph,'Color'))
    end
    function addStatSign(src,event)
        errorbar(xticks,dat,nanstd(gui.plotResultsIndivData)./sqrt(sum(~isnan(gui.plotResultsIndivData))),'Color',get(ph,'Color'))
    end




    function LineSelected(ObjectH, EventData,act_plot)
        if marked_lines(act_plot)>0
            %             for hc = 1:numel(H)
            %                 set(lineH(hc), 'LineWidth', linewidtlineH(hc));
            %             end
            set(lineH(act_plot), 'LineWidth', linewidth(act_plot));
            marked_lines(act_plot) = 0;
        else
            marked_lines(marked_lines>0) = marked_lines(marked_lines>0) + 1;
            marked_lines(act_plot) = 1;
            marked_lines(marked_lines>no_possiblemarks) = 0;
            for ml=find(~marked_lines)
                set(lineH(ml), 'LineWidth', linewidtlineH(ml));
            end
            for ml=find(marked_lines)
                set(lineH(ml), 'LineWidth', linewidtlineH(ml)+2);
            end
        end
    end

    function PointSelected(ObjectH, EventData,act_plot)
        C = get (get(ObjectH,'Parent'), 'CurrentPoint');
        dst = sqrt((C(1,1)-xdataCell{act_plot}).^2+(C(1,2)-ydataCell{act_plot}).^2);
        [~, idx] = min(dst);
        
        if marked_points(idx,act_plot)>0
            delete(markh(idx,act_plot));
            marked_points(idx,act_plot) = 0;
        else
            
            marked_points(marked_points>0) = marked_points(marked_points>0) + 1;
            
            marked_points(idx,act_plot) = 1;
            for mp = 1:size(marked_points,1)
                for mp2 = 1:size(marked_points,2)
                    if ishandle(markh(mp,mp2)) && marked_points(mp,mp2)>no_possiblemarks
                        delete(markh(mp,mp2));
                        marked_points(mp,mp2) = 0;
                    end
                end
            end
            
            markh(idx,act_plot) = plot(xdataCell{act_plot}(idx),ydataCell{act_plot}(idx),'s','Color',linecolor(act_plot,:));
            uistack(markh(idx,act_plot),'bottom')
        end
    end
    function LineDelete(ObjectH, EventData,act_plot)
        if any(marked_points(:,act_plot)>0)
            delete(markh(ishandle(markh(:,act_plot)),act_plot));
%             marked_points(:,act_plot) = [];
        end
    end
end
