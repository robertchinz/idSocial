function signH = idSocial_auxiliaries_plotSignificance(lineObj,p,slevel,tail,mode,test_type,Xctrs,Yctrs)

signH = [];
if nargin <4 || isempty(tail)
    tail = 'both';
end

if nargin <5 || isempty(mode)
    mode = 'plot2d';
end
if nargin <6 || isempty(test_type)
    test_type = '';
end
if nargin<7 || isempty(Xctrs)
    Xctrs = [];
end
if nargin<8 || isempty(Yctrs)
    Yctrs = [];
end

if isnumeric(mode)
    xticks_single = mode;
    mode = 'single_points';
end
% p = ones(2,2,10)*.03;
% ax = gca;
% hdles = findobj(get(ax,'Children'),'Type','line');
% slevel = 0.05;
% mode = 'plot2d';

markersize = 8;
markersize_map = 8;
linewidth = 3;
fontsize = 12;

no_groups = size(p,1);
no_parts = size(p,3);

if strcmpi(mode,'dynamicMapPolar') || strcmpi(mode,'MapPolar')
    for rw = 1:size(p,3)-1
        for cl = 1:size(p,4)-1
          
            if p(1,2,rw,cl)<slevel
                sign_handle = plot(Xctrs(cl,rw),Yctrs(cl,rw),'k*','MarkerSize',markersize_map);
                pval_struct.p = p(1,2,rw,cl);
                pval_struct.test_type = test_type;
                set(sign_handle,'UserData',pval_struct);
                signH = [signH sign_handle];
                
            end
        end
    end
%     for k=1:numel(signH)
%         addlistener(signH(k), 'MarkerSize', 'PostSet', @(hAxes, eventData) setMarkerSize(hAxes, eventData,signH(k),signH(setxor(1:numel(signH),k))));
%     end
    %     linkprop(signH,{'Color','MarkerSize','MarkerFaceColor','LineWidth','MarkerEdgeColor'})
        uistack(signH,'top');
    
    ax = gca;
    info = get(ax,'Userdata');
    info.marker_positions = {Xctrs,Yctrs};
    set(ax,'Userdata',info);
    
else
    
    ax = get(lineObj(1),'Parent');
    ylim = get(ax,'YLim');
    ax_yscale=diff(ylim);
    
    isIdSocialPlot = false(1,size(lineObj,2));
    colororder = NaN(size(lineObj,2),3);
    for k=1:size(lineObj,2)
        ud = get(lineObj(k),'Userdata');
        isIdSocialPlot(k) = isstruct(ud)&isfield(ud,'funcstring');
    end
    lineObj = lineObj(isIdSocialPlot);
    for k=1:size(lineObj,2)
        if isprop(lineObj(k),'Color')
            colororder(k,:) = get(lineObj(k),'Color');
        elseif isprop(lineObj(k),'FaceColor')
            ax_colororder = get(gca,'Colororder');
            fc = get(lineObj(k),'Facecolor');
            if isnumeric(fc)
                colororder(k,:) = fc;
            else
                colororder(k,:) = ax_colororder(k);
            end
        end
    end
    
%     xticks = get(ax,'XTick');
    xticks = get(lineObj(1),'XData');
    for xa = 2:numel(lineObj)
        xticks = intersect(xticks, get(lineObj(xa),'XData'));
    end
    ypos_signmarker = ylim(2)*ones(size(xticks,2),1);
    ax_ud = get(ax,'Userdata');
    if isstruct(ax_ud) && isfield(ax_ud,'ypos_signmarker') && ~isempty(ax_ud.ypos_signmarker)
        
        ypos_signmarker= ax_ud.ypos_signmarker;
    else
        ypos_signmarker=ypos_signmarker+ax_yscale/20;
        if isempty(ax_ud) || (isstruct(ax_ud) && isfield(ax_ud,'ypos_signmarker')) && ~isempty(ax_ud.ypos_signmarker)
            info.ypos_signmarker = ypos_signmarker;
            set(ax,'Userdata',info);
        end
    end
    
    if strcmp(tail,'both')
        group_index = NaN(nchoosek(no_groups,2)+no_groups,2);
        count = 1;
        for gr1 = 1:no_groups
            for gr2 = gr1:no_groups
                %             group_index=[(1:(no_groups/2))' ((no_groups/2+1):no_groups)'];
                group_index(count,:) = [gr1,gr2];
                count = count + 1;
            end
        end
    else
        group_index=vectorsize2indexcombinations([no_groups no_groups]);
    end
    
    if strcmpi(mode,'plot2d')
        if isfield(get(lineObj(1)),'UData')
            ydat=cellfun(@(x,y) x+y,get(lineObj,'YData'),get(lineObj,'UData'),'UniformOutput',false);
        else
            ydat=get(lineObj,'YData');
        end
        if size(ydat,1)==1 && ~iscell(ydat)
            ydat={ydat};
        end
        ydat=ydat(cellfun(@(x) size(x,2)>1,ydat));
        xdat = get(lineObj,'XData');
        
        
        
        if no_groups>1
            for act_idx=1:size(group_index,1)
                
                group1=group_index(act_idx,1);
                group2=group_index(act_idx,2);
                if group1 ~= group2
                    if strcmp(tail,'both') % For 'both' p(a,b) should be the same as p(b,a), but it might only be stored at one of the idx combinations.
                        p_act = max(squeeze(p(group1,group2,:)),squeeze(p(group2,group1,:)));
                    elseif strcmp(tail,'left') || strcmp(tail,'right')
                        p_act = squeeze(p(group1,group2,:));
                    end
                    %                 ypos_signmarker=ypos_signmarker+double(squeeze(p(group1,group2,:))<slevel)*ax_yscale/20;
                    ypos_act=ypos_signmarker;
                    ypos_act(p_act>=slevel | isnan(p_act))=NaN;
                    sign_handle=plot(xticks, ...
                        ypos_act,'o','MarkerEdgeColor',colororder(group1,:),'MarkerFaceColor',colororder(group2,:),'MarkerSize',markersize,'LineWidth',linewidth);
                    signH = [signH sign_handle];
                    % Save handle so marker can be deleted if line is deleted
                    line_ud = get(lineObj(group1),'UserData');
                    if ~isempty(line_ud) && isstruct(line_ud) && isfield(line_ud,'significanceMarkerH')
                        line_ud.significanceMarkerH = [line_ud.significanceMarkerH sign_handle];
                    elseif isempty(line_ud) || (isstruct(line_ud) && ~isfield(line_ud,'significanceMarkerH'))
                        line_ud.significanceMarkerH = sign_handle;
                    end
                    set(lineObj(group1),'UserData',line_ud);
                    
                    line_ud = get(lineObj(group2),'UserData');
                    if ~isempty(line_ud) && isstruct(line_ud) && isfield(line_ud,'significanceMarkerH')
                        line_ud.significanceMarkerH = [line_ud.significanceMarkerH sign_handle];
                    elseif isempty(line_ud) || (isstruct(line_ud) && ~isfield(line_ud,'significanceMarkerH'))
                        line_ud.significanceMarkerH = sign_handle;
                    end
                    set(lineObj(group2),'UserData',line_ud);
                    
                    pval_struct.p = squeeze(p(group1,group2,:));
                    pval_struct.test_type = test_type;
                    set(sign_handle,'UserData',pval_struct);
                    ypos_signmarker(~isnan(ypos_act)) = ypos_signmarker(~isnan(ypos_act))+ax_yscale/20;
                    ylim(2)=max(ylim(2),max(ypos_act)*1.08);
                end
            end
            
        end
    elseif strcmpi(mode,'single_points')
        if no_groups>1
            for act_idx=1:size(group_index,1)
                
                group1=group_index(act_idx,1);
                group2=group_index(act_idx,2);
                if group1 ~= group2
                    if group1 ~= group2
                        minGr = min(find(xticks_single(group1)>=xticks,1,'last'),find(xticks_single(group2)>=xticks,1,'last'));
                        maxGr = max(find(xticks_single(group1)<=xticks,1,'first'),find(xticks_single(group2)<=xticks,1,'first'))-1;
                        
                        
                        %                     maxGr = max(find(xticks_single(group1)==xticks),find(xticks_single(group2)==xticks));
                        if strcmp(tail,'both') % For 'both' p(a,b) should be the same as p(b,a), but it might only be stored at one of the idx combinations.
                            p_act = max(squeeze(p(group1,group2)),squeeze(p(group2,group1)));
                        elseif strcmp(tail,'left') || strcmp(tail,'right')
                            p_act = squeeze(p(group1,group2));
                        end
                        %                     if p_act<0.05
                        
                        yp = max(ypos_signmarker(minGr:maxGr));
                        
                        
                        if p_act>0.05
                            star_string = 'n.s.';
                        elseif p_act <=0.0001
                            star_string = '****';
                        elseif p_act <=0.001
                            star_string = '***';
                        elseif p_act <=0.01
                            star_string = '**';
                        elseif p_act<=0.05
                            star_string = '*';
                        end
                        if ~isnan(p_act)
                            sign_handle = plot([xticks_single(group1) xticks_single(group2)],[yp yp], ...
                                'k','Linewidth',linewidth);
                            text_handle = text((xticks_single(group1)+xticks_single(group2))/2,yp,star_string,'VerticalAlignment','bottom','HorizontalAlignment','center', ...
                                'Fontsize',fontsize);
                            signH = [signH sign_handle text_handle];
                            
                            ypos_signmarker = ones(size(ypos_signmarker))*yp+ax_yscale/20;
                            ylim(2)=max(ylim(2),max(yp)*1.08);
                        end
                        %                     end
                    end
                    
                end
            end
        end
        
    end
    
    
    
    if isempty(ax_ud) || isempty(fieldnames(ax_ud))  || (isstruct(ax_ud) && isfield(ax_ud,'ypos_signmarker'))
        info.ypos_signmarker = ypos_signmarker;
        set(ax,'Userdata',info);
    end
    % set(gca,'YLim',ylim);
    set(gca, 'ylimmode', 'auto')
end
% function setMarkerSize(hAxes, eventData,thisH,markerH)
%     set(markerH,'MarkerSize',get(thisH,'MarkerSize'));
% end
end