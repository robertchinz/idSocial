function [input_data,sphandle, legend_handles]=idSocial_plot(input_data,act_method,plot_mode,plot_options)
sphandle = [];
legend_handles = [];
warning('off','MATLAB:gui:latexsup:UnableToInterpretLaTeXString')
warning('off','MATLAB:gui:latexsup:UnsupportedFont')
%% Check input
if nargin<3 || ~isfield(plot_mode,'visible')|| isempty( plot_mode.visible)
    visible='On';
else
    visible=plot_mode.visible;
end
if (nargin<3 || isempty(plot_mode) ) && isfield(input_data.(act_method),'plot_mode')
    plot_mode = input_data.(act_method).plot_mode;
end
if nargin<4 || isempty(plot_options)
    plot_options = [];
end

if  (~isfield(plot_mode,'display_mode') || isempty( plot_mode.display_mode)) && (~isfield(plot_options,'display_mode') || isempty( plot_options.display_mode))
    display_mode='plot2d';
elseif isfield(plot_options,'display_mode') && ~isempty( plot_options.display_mode)
    display_mode=plot_options.display_mode;
elseif isfield(plot_mode,'display_mode') && ~isempty( plot_mode.display_mode)
    display_mode=plot_mode.display_mode;
end
% plot_mode.display_mode='extern';
if ~(isfield(plot_mode,'display_mode') && strcmpi(plot_mode.display_mode,'extern'))
    %%
    % Set options
    options=input_data(1,1,1).(act_method).options;
    [~,options]=idSocial_readparams([],plot_options,options,act_method,1);
    % act_method=         options.act_method;
    edges=              options.edges;
    titlefontsize=      options.plot_titlefontsize;
    %     legendstring =      options.plot_legendstring;
    legendstring = {''};
    tickfontsize=       options.plot_tickfontsize;
    labelfontsize=      options.plot_labelfontsize;
    fontsize_legend =   labelfontsize;
    axis_scaling=       options.plot_axis_scaling;
    colororder=               options.plot_colororder;
    linewidth=          options.plot_linewidth;
    timeintervals_in_min=options.timeintervals_in_min;
    axislabelstring=    options.plot_axislabelstring;
    if isfield(options,'plot_titlestring')
        titlestring=        options.plot_titlestring;
    else
        titlestring=  {''};
    end
    info=               input_data(1,1).info;
    framerate=          info.framerate;
    plot_xticklabel=        options.plot_xticklabel;
    linestyle = options.plot_linestyle;
    plot_bootstrap_repetitions= options.plot_bootstrap_repetitions;
    bootstrap_statfunc=options.plot_bootstrap_statfunc;
    plot_deviation = options.plot_deviation;
    random_data = options.random_data;
    
    
    
    [~, project_path]= idSocial_recursiveGetOptionsFromOptionsCell(input_data(1,1,1).options,'project_path');
    
    significance_between_groups= options.significance_between_groups;
    statistical_test_type=      options.statistical_test_type;
    statistical_test_significance=                  options.statistical_test_significance;

    group_names=    info.group_name;
    
    %% Some standard parameters
    output_plot=input_data(1,1).(act_method).output_plot;
    
    no_subplots=size(output_plot.data,1);
    no_legends=size(output_plot.data,2);

    if isfield(options,'plot_linestyle')
        linestyle = options.plot_linestyle;
    else
        linestyle = repmat('-',[no_legends ,1]);
    end
    
    linestyle = repmat(linestyle,[10,1]);
    if random_data
            linestyle(1:no_legends/2) = repmat(linestyle(1),[no_legends-no_legends/2 ,1]);
            linestyle(no_legends/2+1:no_legends) = repmat({'--'},[no_legends-no_legends/2 ,1]);
    end
   
    duration=           info.duration;
    if isempty(timeintervals_in_min) || timeintervals_in_min> min(min(floor(min(info.no_frames./info.framerate)/60)));
        timeintervals_in_min=max(duration(duration>0));
    end;
    
    
    no_frames_part_array=floor(timeintervals_in_min.*framerate*60)-1;
    no_parts=2*floor(info.no_frames./no_frames_part_array)+1;
    min_no_parts=min(no_parts(:));
    max_no_parts=max(no_parts(:));
    
    fname=mfilename;
    fname=fname(15:end);
    % timevec=[((1:max_no_parts)-1)*timeintervals_in_min/2 inf];
    timevec=[ [((1:max_no_parts)-1)*timeintervals_in_min/2]' ...
        [((1:max_no_parts)-1)*timeintervals_in_min/2]'+timeintervals_in_min];
    timevec=vertcat([0 timeintervals_in_min/2],timevec);
    
    if iscell(edges)
        edges=edges{1};
    end
    
    
    
    %% Plot the stuff
    
    timeidx=find(strcmpi(output_plot.dim_names_original,'time'));
    groupidx=find(strcmpi(output_plot.dim_names_original,'group'));
    edgesXidx=find(strcmpi(output_plot.dim_names_original,'edgeX'));
    edgesYidx=find(strcmpi(output_plot.dim_names_original,'edgeY'));
    
    
    %% Set xaxis
    
    if output_plot.xaxis_idx==edgesXidx && strcmpi(display_mode,'hist')
        fidx=1;
        de=diff(edges);
        xvec=1:output_plot.xaxis_length;
        xticlabel=edges(1:end);%+de(1)/2;
    elseif output_plot.xaxis_idx==edgesXidx && ~strcmpi(display_mode,'hist')
        fidx=1;
        de=diff(edges);
        xvec=1:output_plot.xaxis_length;
        xticlabel=edges(1:end-1)+de(1)/2;
        
        
    elseif output_plot.xaxis_idx==timeidx
        xvec=1:output_plot.xaxis_length;
        xticlabel=(xvec-1)*timeintervals_in_min/2;
        axislabelstring{1}='Time [min]';
    else
        xvec=1:output_plot.xaxis_length;
        if isempty(plot_xticklabel)
            xticlabel=xvec;
        else
            xticlabel=plot_xticklabel;
        end
        
        if isempty(axislabelstring{1})
            axislabelstring{1}='';
        end
    end
    if numel(xvec)>1
            dx=xvec(2)-xvec(1);
        else
            dx=1;
    end
    % if numel(xvec)==1
    %     xvec=[xvec xvec];
    %     xticlabel=[xticlabel xticlabel];
    % end
    %% Calculate y- (or color-) limits
    lim=[0 0];
    
    
   
    
    v=[];
    v_norm=[];
    v_dev_up=[];
    v_dev_down=[];
    conf_intervall=cell(no_subplots,no_legends,2);
    for sp=1:no_subplots
        for lg=1:no_legends
            if numel(xticlabel)>1
                tick_diff=(xticlabel(2)-xticlabel(1));
            else
                tick_diff=1;
            end
            v=[v output_plot.data{sp,lg}(:)'];
            v_norm=[v_norm output_plot.data{sp,lg}(:)'/(nansum(output_plot.data{sp,lg}(:)) *tick_diff)];
            if ~isempty(output_plot.data_sign) && ~isempty(output_plot.data_sign{sp,lg}) && plot_deviation
                switch output_plot.statistics_type{sp,lg}
                    case {'Mean','mean','MEAN'}
                        v_dev_up=[v_dev_up output_plot.data_dev{sp,lg}(:)'./sqrt(output_plot.data_no_datapoints{sp,lg}(:)')];
                        v_dev_down=v_dev_up;
                    case {'Median','median','MEDIAN'}
                        S.type='()';
                        S.subs=repmat({':'},[size(output_plot.cell_size,2),1]);
                        S.subs{output_plot.statistics_on_idx(sp,lg)}=2;
                        conf_intervall{sp,lg,2}=squeeze(subsref(output_plot.data_dev{sp,lg},S));
                        v_dev_up=[v_dev_up conf_intervall{sp,lg,2}(:)'];
                        
                        S.subs=repmat({':'},[size(output_plot.cell_size,2),1]);
                        S.subs{output_plot.statistics_on_idx(sp,lg)}=1;
                        conf_intervall{sp,lg,1}=subsref(output_plot.data_dev{sp,lg},S);
                        v_dev_down=[v_dev_down conf_intervall{sp,lg,1}(:)'];
                end
            else
                v_dev_down=0;
                v_dev_up=0;
            end
            
        end
        
    end
    if isempty(v_dev_down); v_dev_down=NaN; end;
    if isempty(v_dev_up); v_dev_up=NaN; end;
    maxv=max(v(:)+v_dev_up(:));
    minv=min(v(:)-v_dev_down(:));
    maxv_norm=max(v_norm(:));
    
    switch display_mode
        case 'plot2d'
            
            range_ratio=(max(v(:))-min(v(:)))/min(v(:));
            
            if min(lim(1),min(v(:)))>=0
                lim=[0 max(lim(2),max(v(:)+v_dev_up(:)))*1.1];
            else
                lm=max(max(abs(lim)),max(max([abs(v(:)+v_dev_up(:)),abs(v(:)-v_dev_down(:))])));
                lim=[-max(max(abs(lim)),lm) max(max(abs(lim)),lm)]*1.1 ;
            end
        case 'prob'
            
            
            lim=[minv*.9 max(lim(2),maxv)*1.1];
        case 'hist'
            
            lim=[0 maxv_norm*1.1];
    end
    if lim(1)==lim(2); lim(1)=lim(1)-.1;  lim(2)=lim(2)+.1; end
    
    %%
    figname=act_method;
    titlestr_prefix=[];
    fh_raw=figure('Visible',visible,'units','normalized','Color',[.95,.95,.95],...
        'Position',[0.1 0.1 .8 .84]);
    set(fh_raw,'Tag',[act_method '_raw'])
    set(fh_raw,'Name', figname);
    set(fh_raw, 'renderer', 'zbuffer');
    
    sprows=round(sqrt(no_subplots));
    spcols=ceil(sqrt(no_subplots));
    
    def_height_txt=.2;
    % hp=axes('Units','normalized','OuterPosition',[0 def_height_txt 1 1-def_height_txt]);
    
    annotstr='Figure shows: ';
    subplotstring={''};
    sphandle=NaN(no_subplots,1);
%     legend_handles = struct(no_subplots,1);
    for sp=1:max(1,no_subplots)
        row_act=sprows-floor((sp-1)/spcols+1)+1;
        col_act=mod((sp-1),spcols)+1;
        
        sphandle(sp)=axes('Units','normalized','OuterPosition',[(col_act-1)*1/spcols def_height_txt+(row_act-1)*(1-def_height_txt)/sprows 1/spcols (1-def_height_txt)/sprows]);
        hold on
        if ~isempty(output_plot.subplot_idx)
            sps_act=output_plot.subplot_idx(sp,:);
            
            
            for ss=1:length(output_plot.subplot_idx{sp,1})
                if output_plot.subplot_idx{sp,1}(ss)==timeidx
                    output_plot.subplotstring{sp}=['t= ' num2str(timevec(output_plot.subplot_idx{sp,2}(ss),1),'%.1f') '-' num2str(timevec(min(output_plot.subplot_idx{sp,2}(ss),length(timevec)),2),'%.1f') 'min'];
                elseif output_plot.subplot_idx{sp,1}(ss)==groupidx
                    output_plot.subplotstring{sp}=group_names{sps_act{ss,2}};
                end
            end
        end
        
        no_legends=size(output_plot.data,2);
        %     if no_subplots>1
        %         annotstr=[annotstr '#\textbf{' num2str(sp) '.) }'];
        %     end
        for lg=1:no_legends
            %         lg_act=output_plot.legend_idx{sp};
            hold on
            switch display_mode
                case {'plot2d','prob'}
                    xvec_errorbar=xvec+((lg-1)-.5*(no_legends-1))/(.5*(no_legends))*.05*axis_scaling(1);
                    if (strcmpi(output_plot.statistics_type{sp,lg},'MEAN') || strcmpi(output_plot.statistics_type{sp,lg},'MODE') ) && ~isempty(output_plot.data_sign) ...
                            && ~all(cellfun(@(x) isempty(x),output_plot.data_sign(:))) ...
                            && plot_deviation
                        axes(sphandle(sp))
                        
                        he= errorbar(xvec_errorbar,...
                            squeeze(output_plot.data{sp,lg}),...
                            squeeze(output_plot.data_dev{sp,lg})./squeeze(sqrt(output_plot.data_no_datapoints{sp,lg})),linestyle{lg},...
                            'MarkerSize',6,'Color',colororder(lg,:),'linewidth',linewidth);
                        %                      he= errorbar(xvec_errorbar,...
                        %                         squeeze(output_plot.data{sp,lg}),...
                        %                         squeeze(output_plot.data_dev{sp,lg}),...
                        %                         'MarkerSize',6,'Color',colororder(lg,:),'linewidth',linewidth);
                        set(get(he,'Children'),{'LineWidth'},{3; 1})
                    elseif strcmpi(output_plot.statistics_type{sp,lg},'MEDIAN') && ~isempty(output_plot.data_sign)...
                            && plot_deviation
                        
                        he= errorbar(xvec_errorbar,...
                            squeeze(output_plot.data{sp,lg}),...
                            squeeze(conf_intervall{sp,lg,1}),squeeze(conf_intervall{sp,lg,2}),linestyle,...
                            'MarkerSize',6,'Color',colororder(lg,:),'linewidth',linewidth);
                        set(get(he,'Children'),{'LineWidth'},{3; 1})
                    elseif isempty(output_plot.data_sign) || all(cellfun(@(x) isempty(x),output_plot.data_sign(:))) || ~plot_deviation
                        he= plot(xvec_errorbar,...
                            squeeze(output_plot.data{sp,lg}),linestyle,...
                            'MarkerSize',6,'Color',colororder(lg,:),'linewidth',linewidth);
                    else
                        he= plot(xvec,...
                            squeeze(output_plot.data{sp,lg}),...
                            'MarkerSize',6,'Color',colororder(lg,:),'linewidth',linewidth);
                    end
                    
                case  {'hist';'density'}
                    
                    tick_ratio=1;%max(floor(no_ticks_manual/no_ticks),1);
                    zeropos=find(xticlabel==0);
                    if ~isempty(zeropos)
                        xvec_1sthalf=xvec(zeropos:-tick_ratio:1); xvec_1sthalf=xvec_1sthalf(end:-1:1);
                        xvec_2ndhalf= xvec(zeropos+tick_ratio:tick_ratio:end);
                        xvec_new=[xvec_1sthalf xvec_2ndhalf];
                    else
                        xvec_new=xvec(1:tick_ratio:end);
                    end
                    xticlabel_new=xticlabel(xvec_new);
                    switch display_mode
                        case 'density'
                            
                            he= stairs(xticlabel_new,...
                                squeeze(output_plot.data{sp,lg})/(nansum(output_plot.data{sp,lg})*(xticlabel(2)-xticlabel(1))),...
                                'MarkerSize',6,'Color',colororder(lg,:),'linewidth',linewidth);
                        case 'hist'
                            he= stairs(xticlabel_new,...
                                squeeze(output_plot.data{sp,lg}),...
                                'MarkerSize',6,'Color',colororder(lg,:),'linewidth',linewidth);
                    end
                    set(gca,'FontSize',tickfontsize);
            end
            
        end
        no_ticks=length(get(gca,'XTick'));
        
        % Y-axis at x=0;
        if ischar(xticlabel(1))
            xlim=[str2double(xticlabel(1))  str2double(xticlabel(end-1))];
        else
            if numel(xticlabel)>1
                xlim=[xticlabel(1)  xticlabel(end-1)];
            else
                xlim=[xticlabel(1) xticlabel(1)];
            end
        end
        
        
        xlabel(axislabelstring{1},'fontsize',labelfontsize)
        ylabel(axislabelstring{2},'fontsize',labelfontsize);
        try
            
            if ~isempty(output_plot.subplotstring)
                title([titlestr_prefix titlestring{1} ' ' output_plot.subplotstring{sp}],'fontsize',titlefontsize);
            else
                title([titlestr_prefix titlestring{1} ],'fontsize',titlefontsize);
            end
        catch
            keyboard
        end
        
        if significance_between_groups && ~isempty(output_plot.data_sign) %&& strcmpi(output_plot.statistics_type{sp,lg},'MEAN')
            
            no_data_points=NaN(no_legends,1);
            for legendentr=1:no_legends
                if ~strcmpi(output_plot.statistics_type{legendentr},'Hist')
                    no_data_points(legendentr)=size(output_plot.data_sign{sp,legendentr},output_plot.statistics_on_idx(sp,legendentr));
                else
                    no_data_points(legendentr)=size(output_plot.data_signMedian{sp,legendentr},output_plot.statistics_on_idx(sp,legendentr));
                end
            end
            max_data_points=max(no_data_points);
            if strcmpi(display_mode,'hist')
                signval_all=NaN(no_legends, max_data_points);
            else
                signval_all=NaN(no_legends, max_data_points,output_plot.xaxis_length);
            end
            for legendentr=1:no_legends
                if ~isempty(output_plot.data_sign{sp,legendentr})
                    if strcmpi(display_mode,'hist')
                        signval=output_plot.data_signMedian{sp,legendentr};
                        try
                            signval_all(legendentr,1:no_data_points(legendentr),:)=signval;%/(xticlabel(2)-xticlabel(1))+xticlabel(1);
                        catch
                            keyboard
                        end
                    else
                        signval=output_plot.data_sign{sp,legendentr};
                        sign_dims=ndims(signval);
                        sign_indices=[output_plot.statistics_on_idx(sp,legendentr) output_plot.xaxis_idx setdiff(1:sign_dims,[output_plot.statistics_on_idx(sp,legendentr) output_plot.xaxis_idx])];
                        sign_val_perm=squeeze(permute(signval,sign_indices));
                        
                        signval_all(legendentr,1:no_data_points(legendentr),:)=sign_val_perm;
                    end
                end
            end
            disp('Statistical significance is being calculated...');
            [H, p]=idSocial_auxiliaries_StatSign4groups(signval_all,statistical_test_type,statistical_test_significance,gca,colororder,plot_bootstrap_repetitions,display_mode,bootstrap_statfunc);
            disp('Done.')
            input_data(1,1).(act_method).output_plot.pvalue{sp}=p;
            cursorMode = datacursormode(gcf);
%             set(cursorMode, 'enable','on', 'UpdateFcn',@setDataTipTxt, 'NewDataCursorOnClick',false);
            set(cursorMode, 'enable','on', 'UpdateFcn',@setDataTipTxt);
        end
        
        
        % Annotations and legend
        if max(1,no_subplots)>1
            annotstr=[annotstr 'Fig. ' num2str(sp) ') '];
        end
        for lg=1:no_legends
            annotstr=[annotstr output_plot.data_string{sp,lg} '; '];
        end
        
        if significance_between_groups && ~isempty(output_plot.data_sign)
            if strcmpi(output_plot.statistics_type{sp,lg},'MEAN')
                annotstr=[annotstr 'Statistical test: ' statistical_test_type ' at significance level ' num2str(statistical_test_significance)];
            elseif strcmpi(output_plot.statistics_type{sp,lg},'MEDIAN')
                annotstr=[annotstr 'Confidence interval 95%';];
                
            end
        end
        
        
        
        if ~all(isempty([legendstring{:}]))
            
            if size(legendstring,1)>1 && size(legendstring,2) ==1
                legendstring = legendstring';
            end
            [LEGH,OBJH,OUTH,OUTM]=legend(legendstring(sp,:),'Location','Best','FontSize',fontsize_legend);
            set(findobj(OBJH,'-property','LineWidth'),'LineWidth',linewidth)
            legend_handles(sp).LEGH=LEGH;
            legend_handles(sp).OBJH=OBJH;
            legend_handles(sp).OUTH=OUTH;
            legend_handles(sp).OUTM=OUTM;
        
        elseif ~all(isempty([output_plot.legendstring{:}]))
        
            [LEGH,OBJH,OUTH,OUTM]=legend(output_plot.legendstring(sp,:),'Location','Best','FontSize',fontsize_legend);
            set(findobj(OBJH,'-property','LineWidth'),'LineWidth',linewidth)
            legend_handles(sp).LEGH=LEGH;
            legend_handles(sp).OBJH=OBJH;
            legend_handles(sp).OUTH=OUTH;
            legend_handles(sp).OUTM=OUTM;
        end
        
        
        if xlim(1)<0 && xlim(2)>0 && ~strcmp(display_mode,'hist') && ~strcmp(display_mode,'density')
            new_lim=get(gca,'YLim');
            hyax=plot(abs(xlim(1))/(xlim(2)-xlim(1))*(xvec(end)-xvec(1))+(xvec(2)-xvec(1))/2*[1 1],new_lim,'k-','LineWidth',1);
            %         hyax=plot(abs(xlim(1))/(xlim(2)-xlim(1))*(xvec(end)-xvec(1))+xvec(1)*[1 1],new_lim,'k-','LineWidth',1);
            hsAnno = get(hyax, 'Annotation');
            hsLegend = get(hsAnno, 'LegendInformation');
            set(hsLegend, 'IconDisplayStyle', 'off');
        end
        % Annotations for different display modes:
        % X-axis at y=0;
        switch display_mode
            case 'plot2d'
                
                if lim(1)<0 && lim(2)>0
                    xlim=get(gca,'XLim');
                    hxax=plot(xlim,[0 0],'k-','LineWidth',1);
                    hsAnno = get(hxax, 'Annotation');
                    hsLegend = get(hsAnno, 'LegendInformation');
                    set(hsLegend, 'IconDisplayStyle', 'off');
                end
            case 'prob'
                
                xlim=get(gca,'XLim');
                hxax=plot(xlim,[0.5 0.5],'k-','LineWidth',1);
                hsAnno = get(hxax, 'Annotation');
                hsLegend = get(hsAnno, 'LegendInformation');
                set(hsLegend, 'IconDisplayStyle', 'off');
            case 'hist'
                
        end
        box on;
        
        if ~strcmpi(display_mode,'hist') && ~strcmpi(display_mode,'density')
            set(gca,'XLim',[xvec(1)-dx/2 xvec(end)+dx/2]);
            no_ticks_manual=length(xvec);
            tick_ratio=max(floor(no_ticks_manual/no_ticks),1);
            zeropos=find(xticlabel==0);
            if ~isempty(zeropos) && zeropos>1
                xvec_1sthalf=xvec(zeropos:-tick_ratio:1); xvec_1sthalf=xvec_1sthalf(end:-1:1);
                xvec_2ndhalf= xvec(zeropos+tick_ratio:tick_ratio:end);
                xvec_new=[xvec_1sthalf xvec_2ndhalf];
            else
                % Find tick ratio which fits number of ticks:
                pos_tick_list=(1:2*tick_ratio+1);
                act_no_ticks=(length(xvec)-1)./tick_ratio;
                pos_no_ticks=(length(xvec)-1)./pos_tick_list;
                int_no_ticks_idx=find(pos_no_ticks==round(pos_no_ticks));
                tick_dif=pos_no_ticks-act_no_ticks;
                [~, sort_tick_dif_idx]=sort(abs(tick_dif));
                [~,idx]=ismember(sort_tick_dif_idx,int_no_ticks_idx);
                idx= idx(idx~=0);
                idx=int_no_ticks_idx(idx(1));
                tick_ratio_new=pos_tick_list(idx);
                
                xvec_new=xvec(1:tick_ratio_new:end);
            end
            xticlabel_new=xticlabel(xvec_new);
            set(gca,'XTick',xvec_new,'XTickLabel',xticlabel_new,'FontSize',tickfontsize)
        end
        try
            set(sphandle(sp),'UserData','top')
        catch
            keyboard
        end
        
        % Look for matlab-like properties
        objFNames = fieldnames(get(sphandle(sp)));
        optFNames =  fieldnames(options);
        for k = 1:size(objFNames,1)
            cmpStr = strfind(optFNames,objFNames{k});
            strIdx = cellfun(@(x) ~isempty(x), cmpStr);
            if any(strIdx) 
                try
                set(sphandle(sp),objFNames{k},options.(optFNames{strIdx}));
                catch
                    warning([mfilename ': Failed to set ' objFNames{k}])
                end
            end
        end
        
    end
    
    set(fh_raw,'ResizeFcn',@(one,two) resize_annotations(one,two,annotstr))
    ha=axes('Units','normalized','OuterPosition',[0 0 1 def_height_txt],'XTick',[],'YTick',[],'xticklabel',[],'yticklabel',[]);%,'Visible','off');
    set(ha,'UserData','bottom')
    box on
    
    set(fh_raw,'UserData',annotstr);
    txt=text(0.01,0.02, annotstr  ,'Units','normalized','VerticalAlignment','bottom','HorizontalAlignment','left');
    
    % The following is a hack to call 'resize_annotations'
    fgpos=get(fh_raw,'Position');
    set(fh_raw,'Position',[fgpos(1) fgpos(2) fgpos(3)*.99 fgpos(4)])
    
%     idSocial_saveFigure([act_method '_' fname],project_path)
end
end

function output_txt = setDataTipTxt(~,event_obj)
% ~            Currently not used (empty)
% event_obj    Object containing event data structure
% output_txt   Data cursor text (string or cell array
%              of strings)

pos=get(event_obj,'Position');
pvals=get(get(event_obj,'Target'),'UserData');
if ~isempty(pvals)
    output_txt={['x=' num2str(pos(1))];['y=' num2str(pos(2))];['p=' num2str(pvals(pos(1)))]};
else
    output_txt={['x=' num2str(pos(1))];['y=' num2str(pos(2))];};
    
end
end

function resize_annotations(one, two,str)
% kaka
str=get(one,'UserData');
ax1 = findobj(one,'Type','axes','-and','UserData','top');
ax2 = findobj(one,'Type','axes','-and','UserData','bottom');
if ~isempty(ax1) && ~isempty(ax2)
    txt = findobj(one,'Type','text','-and','Parent',ax2);
    %     set(txt,'String',['\textsf{' str '}'])
    set(txt,'String', str)
    extnt=get(txt,'Extent');
    txtpos=get(txt,'Position');
    extntx=extnt(3);
    
    len=max(1,length(str)-10);%-length('\textsf{}________________');
    spce=strfind(str,' ');
    new_len=floor(len/extntx);
    def_height_txt=get(ax2,'OuterPosition');
    def_height_txt=def_height_txt(4);
    
    break_string=str;
    spce_diff=[spce length(str)];
    for k=1:length(spce)
        if spce_diff(k+1)>new_len
            break_string(spce(k))='#';
            
            spce_diff=spce_diff-spce_diff(k);
        end
    end
    
    break_string=textscan(break_string,'%s','Delimiter','#');
    %     set(txt,'String', strcat('\textsf{' ,break_string{1},'}'));
    set(txt,'String', break_string{1});
    
    %     set(txt,'EdgeColor','red')
    extnt=get(txt,'Extent');%
    extnty=extnt(4);
    set(ax2,'OuterPosition',[0 0 1 def_height_txt*extnty])
    
    for k=1:size(ax1,1)
        act_pos=get(ax1(k),'OuterPosition');
        posy=act_pos(2);
        posx=act_pos(1);
        width=act_pos(3);
        new_height=act_pos(4)/(1-def_height_txt)*(1-def_height_txt*extnty);
        
        new_posy =  def_height_txt*extnty+(posy-def_height_txt)/(1-def_height_txt) * (1-def_height_txt*extnty);
        
        %         new_posy=def_height_txt*extnty+(posy-def_height_txt)*(1-def_height_txt*extnty);%/(1-def_height_txt)*extnty;
        
        set(ax1(k),'OuterPosition',[posx new_posy width new_height])
    end
end
end