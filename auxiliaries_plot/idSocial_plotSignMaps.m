function input_data=idSocial_plotSignMaps(input_data,act_method,plot_mode)

if nargin<3 || ~isfield(plot_mode,'visible')|| isempty( plot_mode.visible)
    visible='On';
else
    visible=plot_mode.visible;
end

load('colormap4InteractionMaps.mat'); % Load colormap
if ~(isfield(plot_mode,'display_mode') && strcmpi(plot_mode.display_mode,'extern'))
    % Set options
    options=input_data(1,1).(act_method).options;
    
    discard_percentage_outliers=options.discard_percentage_outliers; % This will cut off the upper and bottom x% of the z-values when calculating the color limits for the maps.
    
    
    titlefontsize=      options.plot_titlefontsize;
    tickfontsize=       options.plot_tickfontsize;
    axislabelsize=      options.plot_axislabelsize;
    labelfontsize=      options.plot_labelfontsize;
    rowlabelfontsize=   12;%options.rowlabelfontsize;
    columnlabelfontsize=  12;  %options.columnlabelfontsize;
    gr_tr_fontsize=     options.plot_gr_tr_fontsize;
    colororder=               options.plot_colororder;
    spacing=            options.spacing;
    linewidth=          options.plot_linewidth;
    [~, project_path]= idSocial_recursiveGetOptionsFromOptionsCell(input_data(1,1,1).options,'project_path');
    timeintervals_in_min=options.timeintervals_in_min;
    axislabelstring=    options.plot_axislabelstring;
    titlestring=        options.plot_titlestring;
    info=               input_data(1,1).info;
    framerate=          info.framerate;
    significance_between_groups= options.significance_between_groups;
    statistical_test_type=      options.statistical_test_type;
    statistical_test_significance=                  options.statistical_test_significance;
    plot_bootstrap_repetitions=  options.plot_bootstrap_repetitions;
    plot_mode=      options.plot_mode;

    
    if isfield(options,'plot_caption')
        caption   = options.plot_caption;
    else
        caption=true;
    end
    
    if isfield(options,'plot_clims')
        plot_clims   = options.plot_clims;
    else
        plot_clims=[];
    end
    
    %% Some standard parameters
    output_plot=input_data(1,1).(act_method).output_plot;
    
    no_groups=size(input_data,1);
    trials=info.trials;
    max_no_frames=max(info.no_frames);
    duration=           info.duration;
    
    if isempty(timeintervals_in_min) || timeintervals_in_min> min(min(floor(min(info.no_frames./info.framerate))/60));
        timeintervals_in_min=max(max(duration(duration>0)));
    end;
    
    if caption
        def_height_txt=.2;
    else
        def_height_txt = 0;
    end
    no_frames_part_array=floor(timeintervals_in_min.*framerate*60)-1;
    no_parts=2*floor(info.no_frames./no_frames_part_array)+1;
    min_no_parts=min(min(min(no_parts)));
    max_no_parts=max(max(max(no_parts)));
    
    group_names=    info.group_name;
    edgesx=options.edges{1};
    ticsx=floor(edgesx(1)):ceil(edgesx(end));
    edgesy=options.edges{2};
    ticsy=floor(edgesy(1)):ceil(edgesy(end));
    fname=mfilename;
    fname=fname(15:end);
    timevec=[ [((1:max_no_parts)-1)*timeintervals_in_min/2]' ...
        [((1:max_no_parts)-1)*timeintervals_in_min/2]'+timeintervals_in_min];
    timevec=vertcat([0 timeintervals_in_min/2],timevec);
    
    
    
    timeidx=find(strcmpi(output_plot.dim_names_original,'time'));
    groupidx=find(strcmpi(output_plot.dim_names_original,'group'));
    edgesXidx=find(strcmpi(output_plot.dim_names_original,'edgeX'));
    edgesYidx=find(strcmpi(output_plot.dim_names_original,'edgeY'));
    %% Calculate statistical significance
    no_rows=cellfun(@(x) size(x,7),output_plot.data);
    no_cols=cellfun(@(x) size(x,8),output_plot.data);
    
    if ~all(all(no_rows==no_rows(1))) || ~all(all(no_cols==no_cols(1)))
        error([mfilename ': All maps must have the same size!'])
    end
    no_rows=no_rows(1);
    no_cols=no_cols(1);
    
    no_figs=size(output_plot.data,1);
    no_comparemaps=size(output_plot.data,2);
    maps_grid_all_figures=NaN(no_comparemaps+1,no_comparemaps+1,no_rows,no_cols,no_figs);
    H=NaN(no_figs,no_comparemaps+1,no_comparemaps+1,no_rows,no_cols);
    pVal=NaN(no_figs,no_comparemaps+1,no_comparemaps+1,no_rows,no_cols);
    for fg=1:no_figs
        
        for cmpr1=2:no_comparemaps+1
            
            map1=squeeze(output_plot.data{fg,cmpr1-1});
            statmaps1=output_plot.data_sign{fg,cmpr1-1};
            sign_idx=1:ndims(output_plot.data_sign{fg,cmpr1-1});
            sign_idx=[output_plot.statistics_on_idx(fg,cmpr1-1) setxor(output_plot.statistics_on_idx(fg,cmpr1-1),sign_idx)];
            statmaps1=squeeze(permute(statmaps1,sign_idx));
            for cmpr2=cmpr1+1:no_comparemaps+1
                map2=squeeze(output_plot.data{fg,cmpr2-1});
                statmaps2=output_plot.data_sign{fg,cmpr2-1};
                sign_idx=1:ndims(output_plot.data_sign{fg,cmpr2-1});
                sign_idx=[output_plot.statistics_on_idx(fg,cmpr2-1) setxor(output_plot.statistics_on_idx(fg,cmpr2-1),sign_idx)];
                statmaps2=squeeze(permute(statmaps2,sign_idx));
                
                
                maps_grid_all_figures(cmpr1,cmpr2,:,:,fg)= map1-map2;
                if all(arrayfun(@(k) sum(~isnan(statmaps1(:,k))),1:no_rows*no_cols)>1) && ...
                        all(arrayfun(@(k) sum(~isnan(statmaps2(:,k))),1:no_rows*no_cols)>1)
                    switch statistical_test_type
                        case 'bootstrap'
                            [Htemp ptemp]=arrayfun(@(k)idSocial_bootstrap(statmaps1(:,k),statmaps2(:,k),plot_bootstrap_repetitions),1:no_rows*no_cols);
                        case 'ttest'
                            [Htemp ptemp]=arrayfun(@(k)ttest2(statmaps1(:,k),statmaps2(:,k),statistical_test_significance,[],'unequal'),1:no_rows*no_cols);
                        otherwise
                            try
                                [ptemp Htemp]=arrayfun(@(k)ranksum(statmaps1(~isnan(statmaps1(:,k)),k),statmaps2(~isnan(statmaps2(:,k)),k),'alpha',statistical_test_significance),1:no_rows*no_cols);
                            catch
                                keyboard
                            end
                    end
                    
                    H(fg,cmpr1,cmpr2,:,:)=reshape(Htemp,no_rows,no_cols);
                    pVal(fg,cmpr1,cmpr2,:,:)=reshape(ptemp,no_rows,no_cols);
                    H(fg,cmpr2,cmpr1,:,:)=H(fg,cmpr1,cmpr2,:,:);
                    pVal(fg,cmpr2,cmpr1,:,:)=pVal(fg,cmpr1,cmpr2,:,:);
                    maps_grid_all_figures(cmpr2,cmpr1,:,:,fg)=pVal(fg,cmpr2,cmpr1,:,:);%-maps_grid_all_figures(cmpr1,cmpr2,:,:);
                end
            end
            maps_grid_all_figures(1,cmpr1,:,:,fg)=output_plot.data{fg,cmpr1-1};
            maps_grid_all_figures(cmpr1,1,:,:,fg)=output_plot.data{fg,cmpr1-1};
        end
        maps_grid_all_figures(1,1,:,:,fg)=maps_grid_all_figures(1,2,:,:,fg);
    end
    %% Calculate color limits
    rows_subplot=no_comparemaps+1;
    cols_subplot=no_comparemaps+1;
    
    allzvals=[];
    allzvals_norm=[];
    for fg=1:no_figs
        for k=1:no_comparemaps
            allzvals=[allzvals output_plot.data{fg,k}(:)'];
            allzvals_norm=[allzvals_norm (output_plot.data{fg,k}(:)/nansum(output_plot.data{fg,k}(:))/(edgesx(2)-edgesx(1))./(edgesy(2)-edgesy(1)))'];
        end
    end
    
    clim=-inf;
    clim_meanmaps= max(clim,max(abs(prctile(allzvals,discard_percentage_outliers)),abs(prctile(allzvals,100-discard_percentage_outliers))));
    clim_meanmaps_bottom= max(clim,min(abs(prctile(allzvals,discard_percentage_outliers)),abs(prctile(allzvals,100-discard_percentage_outliers))));
    
    clim_norm=-inf;
    clim_meanmaps_norm= max(clim,max(abs(prctile(allzvals_norm,discard_percentage_outliers)),abs(prctile(allzvals_norm,100-discard_percentage_outliers))));
    clim_meanmaps_norm_bottom= max(clim,min(abs(prctile(allzvals_norm,discard_percentage_outliers)),abs(prctile(allzvals_norm,100-discard_percentage_outliers))));
    
    mapsdiff_all=[];
    for fg=1:no_figs
        mapsdiff=permute(maps_grid_all_figures(2:end,2:end,:,:,fg),[3 4 1 2]);
        mapsdiff=mapsdiff(:,:,logical(triu(ones(no_comparemaps,no_comparemaps))));
        mapsdiff_all=[mapsdiff_all mapsdiff(:)'];
    end
    clim= max(clim,max(abs(prctile(mapsdiff_all,discard_percentage_outliers)),abs(prctile(mapsdiff_all,100-discard_percentage_outliers))));
    clim_diffmaps=clim;
    
    
    clims_all=NaN(rows_subplot,cols_subplot,2);
    clims_accum=zeros(rows_subplot*cols_subplot+1,2);
    
    %% Plot
    gl = opengl('data');
    if strfind(gl.Vendor,'ATI')
        opengl('software')
    end
    
    
    edgesx=options.edges{1};
    ticsx=floor(edgesx(1)):ceil(edgesx(end));
    edgesy=options.edges{2};
    ticsy=floor(edgesy(1)):ceil(edgesy(end));
    
    %% Complete subplot/figure titles and legends
    if ~isempty(output_plot.subplot_idx)
        
        for sp=1:no_figs
            sps_act=output_plot.subplot_idx(sp,:);
            for ss=1:length(output_plot.subplot_idx{sp,1})
                if output_plot.subplot_idx{sp,1}(ss)==timeidx
                    output_plot.subplotstring{sp}=['t= ' num2str(timevec(output_plot.subplot_idx{sp,2}(ss),1),'%.1f') '-' num2str(timevec(min(output_plot.subplot_idx{sp,2}(ss),length(timevec)),2),'%.1f') 'min'];
                elseif output_plot.subplot_idx{sp,1}(ss)==groupidx
                    output_plot.subplotstring{sp}=group_names{sps_act{ss,2}};
                end
            end
        end
    end
    
    %%
    for fg=1:no_figs
        annotstr='Figure shows: ';
        if ~isempty(output_plot.subplotstring) && ~isempty(output_plot.subplotstring{fg})
            figname= [act_method ', ' output_plot.subplotstring{fg}];
        else
            figname=act_method;
        end
        
        %     fh_raw=figure('Visible',visible,'units','normalized','outerposition',[.2 .1 .7 .8],'Color',[.95,.95,.95]);
        fh_raw=figure('Visible',visible,'units','normalized','outerposition',[.2 .1 .7 .8],'Color','w');
        
        set(fh_raw,'Tag',[act_method '_raw'])
        
        set(fh_raw,'Name',figname);
        
        
        
        
        pos_subplot=0;
        
        cmap1 = jet(256);
        cmap2 = colormap4InteractionMaps;
        cmap_all=cat(1,cmap1,cmap2);
        
        
        
        % Define a colormap that consists of several separate colormaps.
        % Idea from http://www.mathworks.es/support/solutions/en/data/1-GNRWEH/index.html
        cmap_subplot=NaN(rows_subplot,cols_subplot);
        for cmpr1=1:rows_subplot
            for cmpr2=1:cols_subplot
                pos_subplot=pos_subplot+1;
                map_act=squeeze(maps_grid_all_figures(cmpr1,cmpr2,:,:,fg));
                
                
                if cmpr1>1 && cmpr2>1 && cmpr2>cmpr1 % Differences
                    cmap_subplot(cmpr1,cmpr2)=2;
                    clim=clim_diffmaps;%max(abs(prctile(map_act(:),discard_percentage_outliers)),abs(prctile(map_act(:),100-discard_percentage_outliers)));
                    
                    clims_all(cmpr1,cmpr2,:)=[-clim clim];
                    clims_accum(pos_subplot,:)=[clim 2*clim];
                    
                elseif cmpr1>1 && cmpr2>1 && cmpr2<cmpr1 % P-Value
                    cmap_subplot(cmpr1,cmpr2)=1;
                    clim=.1;%max(map_act(:));
                    clims_all(cmpr1,cmpr2,:)=[0 clim];
                    clims_accum(pos_subplot,:)=[0 clim];
                else % Original maps
                    maxv=max(map_act(:));
                    minv=min(map_act(:));
                    
                    range_ratio=(maxv-minv)/minv;
                    if min(map_act(:))>=0 && range_ratio>.7
                        cmap_subplot(cmpr1,cmpr2)=1;
                        map_act = map_act./sum(map_act(:))./(edgesx(2)-edgesx(1))./(edgesy(2)-edgesy(1));
                        maps_grid_all_figures(cmpr1,cmpr2,:,:,fg)=map_act;
                        clim=clim_meanmaps_norm;%max(map_act(:));
                        if clim<1 && clim>.5
                            clim=1;
                        end
                        clims_all(cmpr1,cmpr2,:)=[0 clim];
                        clims_accum(pos_subplot,:)=[0 clim];
                    elseif min(map_act(:))>=0 && range_ratio<.7 && ~( clim_meanmaps_bottom<.5 && clim_meanmaps>.5)
                        cmap_subplot(cmpr1,cmpr2)=1;
                        %                     clim=clim_meanmaps;%max(map_act(:));
                        clims_all(cmpr1,cmpr2,:)=[clim_meanmaps_bottom clim_meanmaps];
                        clims_accum(pos_subplot,:)=[clim_meanmaps_bottom clim_meanmaps];
                    elseif min(map_act(:))>=0 && range_ratio<.7 && ( clim_meanmaps_bottom<.5 && clim_meanmaps>.5)
                        cmap_subplot(cmpr1,cmpr2)=1;
                        %                     clim=clim_meanmaps;%max(map_act(:));
                        clims_all(cmpr1,cmpr2,:)=[.5-max(abs(.5-clim_meanmaps_bottom),abs(.5-clim_meanmaps)) .5+max(abs(.5-clim_meanmaps_bottom),abs(.5-clim_meanmaps))];
                        clims_accum(pos_subplot,:)=[.5-max(abs(.5-clim_meanmaps_bottom),abs(.5-clim_meanmaps)) .5+max(abs(.5-clim_meanmaps_bottom),abs(.5-clim_meanmaps))];
                    elseif min(map_act(:))<0
                        cmap_subplot(cmpr1,cmpr2)=2;
                        clim=clim_meanmaps;%max(map_act(:));
                        clims_all(cmpr1,cmpr2,:)=[-clim clim];
                        clims_accum(pos_subplot,:)=[clim 2*clim];
                    else
                        cmap_subplot(cmpr1,cmpr2)=1;
                        clim=1;%max(map_act(:));
                        clims_all(cmpr1,cmpr2,:)=[-clim clim];
                        clims_accum(pos_subplot,:)=[clim 2*clim];
                    end
                    
                end
            end
        end
        
        
        
        
        if no_comparemaps > 1
            pos_subplot=(no_comparemaps+1)*(no_comparemaps+1)+1;%cols_subplot;
            subplot_handle=NaN((no_comparemaps+1)*(no_comparemaps+1),1);
            for cmpr1=no_comparemaps+1:-1:1
                for cmpr2=no_comparemaps+1:-1:1
                    pos_subplot=pos_subplot-1;
                    row_act=no_comparemaps+1-cmpr1+1;
                    col_act=cmpr2;
                    if ~isempty(maps_grid_all_figures(cmpr1,cmpr2,:,:,fg))
                        try
                            subplot_handle(pos_subplot)=...
                                axes('Units','normalized','OuterPosition',[(col_act-1)*1/cols_subplot def_height_txt+(row_act-1)*(1-def_height_txt)/rows_subplot 1/cols_subplot (1-def_height_txt)/rows_subplot]);
                            set(gca,'UserData','top')
                        catch
                            keyboard
                        end
                        
                        %                     p = get(subplot_handle(pos_subplot), 'pos');
                        %                     p(2)=p(2) - 0.05;
                        %                     p(3:4) = p(3:4) + 0.02;
                        %                     set(subplot_handle(pos_subplot), 'pos', p);
                        
                        map_act=squeeze(maps_grid_all_figures(cmpr1,cmpr2,:,:,fg));
                        H_act=squeeze(H(fg,cmpr1,cmpr2,:,:));
                        if cmpr1~=cmpr2 && cmpr2>1 && cmpr1>1
                            
                            if no_cols==1;
                                map_act=map_act';
                                H_act=H_act';
                            end
                            %%%%%%%
                            map_act(map_act> clims_all(cmpr1,cmpr2,2))=clims_all(cmpr1,cmpr2,2);
                            map_act(map_act< clims_all(cmpr1,cmpr2,1))=clims_all(cmpr1,cmpr2,1);
                            
                            him= imagesc(edgesx,edgesy,((map_act-clims_all(cmpr1,cmpr2,1))/(clims_all(cmpr1,cmpr2,2)-clims_all(cmpr1,cmpr2,1)))*255+(cmap_subplot(cmpr1,cmpr2)-1)*256);
                            set(gca,'UserData','top')
                            no_xticks=length(get(gca,'XTick')); no_yticks=length(get(gca,'YTick'));
                            colormap(cmap_all)
                            hact=gca;
                            set(gca,'YDir','normal','YtickLabelMode','Manual','YtickMode','Manual')
                            set(gca,'CLim', [0 max(cmap_subplot(:))*256])
                            alphaval=.4;
                            if cmpr2>cmpr1
                                Halpha=double(H_act); Halpha(Halpha==0)=alphaval;
                                set(him,'AlphaDataMapping','none','AlphaData',Halpha)
                            end
                            axis square tight
                            no_xticks_manual=length(ticsx);
                            no_yticks_manual=length(ticsy);
                            
                            if no_xticks_manual==no_yticks_manual
                                no_xticks=max(no_xticks,no_yticks);
                                no_yticks=no_xticks;
                            end
                            tick_ratio=max(1,floor(no_xticks_manual/no_xticks));
                            zeropos=find(ticsx==0);
                            if ~isempty(zeropos)
                                ticsx_1sthalf=ticsx(zeropos:-tick_ratio:1); ticsx_1sthalf=ticsx_1sthalf(end:-1:1);
                                ticsx_2ndhalf= ticsx(zeropos+tick_ratio:tick_ratio:end);
                                ticsx_new=[ticsx_1sthalf ticsx_2ndhalf];
                            else
                                ticsx_new=ticsx(1:tick_ratio:end);
                            end
                            %                     xticlabel_new=xticlabel(ticsx_new);
                            
                            
                            tick_ratio=max(1,floor(no_yticks_manual/no_yticks));
                            zeropos=find(ticsy==0);
                            if ~isempty(zeropos)
                                ticsy_1sthalf=ticsy(zeropos:-tick_ratio:1); ticsy_1sthalf=ticsy_1sthalf(end:-1:1);
                                ticsy_2ndhalf= ticsy(zeropos+tick_ratio:tick_ratio:end);
                                ticsy_new=[ticsy_1sthalf ticsy_2ndhalf];
                            else
                                ticsy_new=ticsy(1:tick_ratio:end);
                            end
                            
                            set(gca,'Xtick',ticsx_new,'Ytick',ticsy_new);
                            set(gca,'XtickLabel',[],'YtickLabel',[])
                            hold on
                            
                            xc=edgesx(1):(edgesx(end)-edgesx(1))/(no_cols-1):edgesx(end);
                            yc=edgesy(1):(edgesy(end)-edgesy(1))/(no_rows-1):edgesy(end);
                            if cmpr2>cmpr1
                                if all(size(Halpha)>1)
                                    
                                    contr=contourc(xc,yc, H_act,[.25 .25]);
                                    contn=1;
                                    while contn+2<length(contr)
                                        if contn==1
                                            plot(contr(1,contn+1:contn+contr(2,contn)),contr(2,contn+1:contn+contr(2,contn)),'Color',[.5 .5 .5],'LineWidth',2);
                                            contn=contn+contr(2,contn)+1;
                                        else
                                            plot(contr(1,contn+1:contn+contr(2,contn)),contr(2,contn+1:contn+contr(2,contn)),'Color',[.5 .5 .5],'LineWidth',2);
                                            contn=contn+contr(2,contn)+1;
                                        end
                                    end
                                end
                            end
                        elseif cmpr1==1 && cmpr2==1
                            
                            % The following plot will be invsible, left the whole stuff in case I want to put something else
                            %                         if no_cols==1;
                            %                             map_act=map_act';
                            %                         end
                            %                         imagesc(edgesx,edgesy,((map_act-clims_all(cmpr1,cmpr2,1))/(clims_all(cmpr1,cmpr2,2)-clims_all(cmpr1,cmpr2,1)))*255+(cmap_subplot(cmpr1,cmpr2)-1)*256);
                            set(gca,'Visible','off')
                            set(gca,'UserData','top')
                            %                         colormap(cmap_all)
                            %                         set(h,'Visible','off')
                            %
                            %
                            %                         set(gca,'Visible','off')
                            
                        elseif cmpr2==1 && cmpr1>1
                            if no_cols==1;
                                map_act=map_act';
                                
                            end
                            map_act(map_act> clims_all(cmpr1,cmpr2,2))=clims_all(cmpr1,cmpr2,2);
                            map_act(map_act< clims_all(cmpr1,cmpr2,1))=clims_all(cmpr1,cmpr2,1);
                            imagesc(edgesx,edgesy,((map_act-clims_all(cmpr1,cmpr2,1))/(clims_all(cmpr1,cmpr2,2)-clims_all(cmpr1,cmpr2,1)))*255+(cmap_subplot(cmpr1,cmpr2)-1)*256);
                            set(gca,'UserData','top')
                            no_xticks=length(get(gca,'XTick')); no_yticks=length(get(gca,'YTick'));
                            colormap(cmap_all)
                            set(gca,'CLim', [0 max(cmap_subplot(:))*256])
                            
                            hact=gca;
                            set(gca,'YDir','normal','YtickLabelMode','Manual','YtickMode','Manual')
                            
                            axis square tight
                            set(gca,'XtickLabel',[],'YtickLabel',[])
                            
                            no_xticks_manual=length(ticsx);
                            no_yticks_manual=length(ticsy);
                            
                            if no_xticks_manual==no_yticks_manual
                                no_xticks=max(no_xticks,no_yticks);
                                no_yticks=no_xticks;
                            end
                            tick_ratio=max(1,floor(no_xticks_manual/no_xticks));
                            zeropos=find(ticsx==0);
                            if ~isempty(zeropos)
                                ticsx_1sthalf=ticsx(zeropos:-tick_ratio:1); ticsx_1sthalf=ticsx_1sthalf(end:-1:1);
                                ticsx_2ndhalf= ticsx(zeropos+tick_ratio:tick_ratio:end);
                                ticsx_new=[ticsx_1sthalf ticsx_2ndhalf];
                            else
                                ticsx_new=ticsx(1:tick_ratio:end);
                            end
                            
                            
                            
                            tick_ratio=max(1,floor(no_yticks_manual/no_yticks));
                            zeropos=find(ticsy==0);
                            if ~isempty(zeropos)
                                ticsy_1sthalf=ticsy(zeropos:-tick_ratio:1); ticsy_1sthalf=ticsy_1sthalf(end:-1:1);
                                ticsy_2ndhalf= ticsy(zeropos+tick_ratio:tick_ratio:end);
                                ticsy_new=[ticsy_1sthalf ticsy_2ndhalf];
                            else
                                ticsy_new=ticsy(1:tick_ratio:end);
                            end
                            
                            
                            set(gca,'Xtick',ticsx_new,'Ytick',ticsy_new,'YtickLabel',ticsy_new,'FontSize',labelfontsize)
                            
                            
                            %                         B = axes('Position',get(hact,'Position'),'YLim',get(hact,'YLim'),'XLim',get(hact,'XLim'));
                            %                         set(gca,'UserData','top')
                            %                         set(B,'Xtick',ticsx_new,'Ytick',ticsy_new,'YtickLabel',[],'XtickLabel',[],'FontSize',labelfontsize)
                            %                         axis tight square
                            %                         uistack(gca,'bottom')
                        elseif cmpr2>1 && cmpr1==1
                            if no_cols==1;
                                map_act=map_act';
                                
                            end
                            map_act(map_act> clims_all(cmpr1,cmpr2,2))=clims_all(cmpr1,cmpr2,2);
                            map_act(map_act< clims_all(cmpr1,cmpr2,1))=clims_all(cmpr1,cmpr2,1);
                            imagesc(edgesx,edgesy,((map_act-clims_all(cmpr1,cmpr2,1))/(clims_all(cmpr1,cmpr2,2)-clims_all(cmpr1,cmpr2,1)))*255+(cmap_subplot(cmpr1,cmpr2)-1)*256);
                            set(gca,'UserData','top')
                            no_xticks=length(get(gca,'XTick')); no_yticks=length(get(gca,'YTick'));
                            set(gca,'CLim', [0 max(cmap_subplot(:))*256])
                            colormap(cmap_all)
                            hact=gca;
                            set(gca,'YDir','normal','YtickLabelMode','Manual','YtickMode','Manual')
                            axis square tight
                            set(gca,'XtickLabel',[],'YtickLabel',[])
                            
                            no_xticks_manual=length(ticsx);
                            no_yticks_manual=length(ticsy);
                            
                            if no_xticks_manual==no_yticks_manual
                                no_xticks=max(no_xticks,no_yticks);
                                no_yticks=no_xticks;
                            end
                            
                            tick_ratio=max(1,floor(no_xticks_manual/no_xticks));
                            zeropos=find(ticsx==0);
                            if ~isempty(zeropos)
                                ticsx_1sthalf=ticsx(zeropos:-tick_ratio:1); ticsx_1sthalf=ticsx_1sthalf(end:-1:1);
                                ticsx_2ndhalf= ticsx(zeropos+tick_ratio:tick_ratio:end);
                                ticsx_new=[ticsx_1sthalf ticsx_2ndhalf];
                            else
                                ticsx_new=ticsx(1:tick_ratio:end);
                            end
                            %                     xticlabel_new=xticlabel(ticsx_new);
                            
                            tick_ratio=max(floor(no_yticks_manual/no_yticks),1);
                            zeropos=find(ticsy==0);
                            if ~isempty(zeropos)
                                ticsy_1sthalf=ticsy(zeropos:-tick_ratio:1); ticsy_1sthalf=ticsy_1sthalf(end:-1:1);
                                ticsy_2ndhalf= ticsy(zeropos+tick_ratio:tick_ratio:end);
                                ticsy_new=[ticsy_1sthalf ticsy_2ndhalf];
                            else
                                ticsy_new=ticsy(1:tick_ratio:end);
                            end
                            %                     yticlabel_new=yticlabel(ticsy_new);
                            
                            
                            set(gca,'Xtick',ticsx_new,'Ytick',ticsy_new,'YtickLabel',[],'FontSize',labelfontsize)
                            
                            
                            %                         B = axes('Position',get(hact,'Position'),'YLim',get(hact,'YLim'),'XLim',get(hact,'XLim'));
                            %                         set(B,'Xtick',ticsx_new,'Ytick',ticsy_new,'YtickLabel',[],'XtickLabel',[],'FontSize',labelfontsize)
                            %                         set(gca,'UserData','top')
                            %                         axis tight square
                            %                         uistack(gca,'bottom')
                        else % Placeholder
                            h=imagesc(edgesx,edgesy,((map_act-clims_all(1,2,1))/(clims_all(1,2,2)-clims_all(1,2,1)))*255+(cmap_subplot(1,2)-1)*256);
                            no_xticks=length(get(gca,'XTick')); no_yticks=length(get(gca,'YTick'));
                            set(gca,'UserData','top')
                            no_xticks_manual=length(ticsx);
                            no_yticks_manual=length(ticsy);
                            
                            if no_xticks_manual==no_yticks_manual
                                no_xticks=max(no_xticks,no_yticks);
                                no_yticks=no_xticks;
                            end
                            tick_ratio=max(floor(no_xticks_manual/no_xticks),1);
                            zeropos=find(ticsx==0);
                            if ~isempty(zeropos)
                                ticsx_1sthalf=ticsx(zeropos:-tick_ratio:1); ticsx_1sthalf=ticsx_1sthalf(end:-1:1);
                                ticsx_2ndhalf= ticsx(zeropos+tick_ratio:tick_ratio:end);
                                ticsx_new=[ticsx_1sthalf ticsx_2ndhalf];
                            else
                                ticsx_new=ticsx(1:tick_ratio:end);
                            end
                            %                     xticlabel_new=xticlabel(ticsx_new);
                            
                            
                            tick_ratio=max(floor(no_yticks_manual/no_yticks),1);
                            zeropos=find(ticsy==0);
                            if ~isempty(zeropos)
                                ticsy_1sthalf=ticsy(zeropos:-tick_ratio:1); ticsy_1sthalf=ticsy_1sthalf(end:-1:1);
                                ticsy_2ndhalf= ticsy(zeropos+tick_ratio:tick_ratio:end);
                                ticsy_new=[ticsy_1sthalf ticsy_2ndhalf];
                            else
                                ticsy_new=ticsy(1:tick_ratio:end);
                            end
                            %                     yticlabel_new=yticlabel(ticsy_new);
                            set(gca,'CLim', [0 max(cmap_subplot(:))*256])
                            colormap(cmap_all)
                            set(h,'Visible','off')
                            set(gca,'Visible','off')
                            set(get(gca,'Title'),'Visible','on')
                            set(get(gca,'YLabel'),'Visible','On')
                            
                            
                        end
                        
                        
                        %%% Put title, labels etc. depending on rows and columns
                        
                        %%%% Titles 1st row and 1st column, 2nd row
                        if cmpr1==1 && cmpr2>1 %|| cmpr1==2 && cmpr2==2
                            
                            title(output_plot.legendstring{fg,cmpr2-1},'fontsize',labelfontsize)
                        elseif cmpr1==2 && cmpr2==1
                            ht=title('Maps','fontsize',labelfontsize);
                            
                        end
                        
                        %%%% Y-ticks for first column
                        if cmpr2==1 || cmpr2==2 && cmpr1==1
                            hylabel=ylabel( options.plot_axislabelstring{2},'fontsize',labelfontsize);
                            set(gca,'Xtick',ticsx_new,'Ytick',ticsy_new,'YtickLabel',ticsy_new,'FontSize',tickfontsize)
                            
                        end
                        
                        %%%% Y-axis label  for last column, first row
                        
                        if cmpr1==1 && cmpr2==cols_subplot && cmpr1<rows_subplot
                            
                            ylabel(gca,'Maps','fontsize',labelfontsize)
                            set(gca,'yaxislocation','right','Ytick',[]);
                            
                            if cmpr2==cols_subplot-1 && cols_subplot==2 % If N=2 put y-label on both sides
                                B = axes('Position',get(hact,'Position'),'YLim',get(hact,'YLim'),'XLim',get(hact,'XLim'));
                                set(B,'yaxislocation','left','Xtick',ticsx_new,'Ytick',ticsy_new,'YtickLabel',ticsy_new,'FontSize',tickfontsize)
                                set(gca,'UserData','top')
                                axis tight square
                                ylabel(B, options.plot_axislabelstring{2},'fontsize',labelfontsize)
                                uistack(gca,'bottom')
                            end
                            
                        end
                        
                        %%%% Y-axis label  for last column
                        
                        if cmpr1>1 && cmpr2==cols_subplot && cmpr1<rows_subplot || cmpr1>1 && cmpr2==cols_subplot-1 && cmpr1==rows_subplot
                            
                            ylabel(gca,output_plot.legendstring{fg,cmpr1-1},'fontsize',labelfontsize)
                            set(gca,'yaxislocation','right','Ytick',[]);
                            
                            
                            if cmpr2==cols_subplot-1 && cols_subplot==2 % If N=2 put y-label on both sides
                                B = axes('Position',get(hact,'Position'),'YLim',get(hact,'YLim'),'XLim',get(hact,'XLim'));
                                set(B,'yaxislocation','left','Xtick',ticsx_new,'Ytick',ticsy_new,'YtickLabel',ticsy_new,'FontSize',tickfontsize)
                                set(gca,'UserData','top')
                                axis tight square
                                ylabel(B, options.plot_axislabelstring{2},'fontsize',labelfontsize)
                                uistack(gca,'bottom')
                            end
                            
                        end
                        
                        %%%% X-axis label last row and last subplot in second last row
                        if cmpr1==rows_subplot || cmpr1==rows_subplot-1 && cmpr2==cols_subplot
                            xlabel(options.plot_axislabelstring{1},'fontsize',labelfontsize)
                            set(gca,'Xtick',ticsx_new,'Ytick',ticsy_new,'XtickLabel',ticsx_new,'FontSize',tickfontsize)
                            
                        end
                        
                        if cmpr1==1 && rows_subplot==2 && cmpr2==2 % If N=2 put y-label on both sides
                            B = axes('Position',get(hact,'Position'),'YLim',get(hact,'YLim'),'XLim',get(hact,'XLim'));
                            set(B,'Xtick',ticsx_new,'Ytick',ticsy_new,'XtickLabel',ticsx_new,'FontSize',tickfontsize)
                            set(gca,'UserData','top')
                            axis tight square
                            xlabel(B, options.plot_axislabelstring{1},'fontsize',labelfontsize)
                            uistack(gca,'bottom')
                            
                        end
                        
                        
                        
                        
                        
                    end
                    %%%% Colorbars
                    
                    if cmpr1==1 && cmpr2==1
                        subplot_handle(pos_subplot)=...
                            axes('Units','normalized','OuterPosition',[(col_act-1)*1/cols_subplot def_height_txt+(row_act-1)*(1-def_height_txt)/rows_subplot 1/cols_subplot (1-def_height_txt)/rows_subplot]);
                        
                        map_act=squeeze(maps_grid_all_figures(1,2,:,:,fg));
                        h=imagesc(edgesx,edgesy,((map_act-clims_all(1,2,1))/(clims_all(1,2,2)-clims_all(1,2,1)))*255+(cmap_subplot(1,2)-1)*256+1);
                        set(gca,'CLim', [0 max(cmap_subplot(:))*256])
                        set(h,'Visible','off')
                        set(gca,'Visible','off')
                        set(gca,'UserData','top')
                        colormap(cmap_all)
                        
                        hc=colorbar;
                        set(hc,'location','southoutside');
                        set(gca,'UserData','top')
                        labelstring=cellstr(num2str((clims_all(1,2,1):diff(clims_all(1,2,:))/4:clims_all(1,2,2))', '%.2g'));
                        set(hc,'XLim', [(cmap_subplot(1,2)-1)*256+1 (cmap_subplot(1,2))*256],...
                            'Xtick',(cmap_subplot(1,2)-1)*256+1:(256-1)/4 :(cmap_subplot(1,2))*256+1,...
                            'XtickLabel', labelstring,'FontSize',tickfontsize);
                        %                 xticklabel_rotate(XTick,rot,XTickLabel
                        title(hc,{titlestring{1}; ' '; 'Means'},'fontsize',labelfontsize)
                        hc_pos=get(hc,'Position');
                        hc_pos(1)=hc_pos(1)+.01;
                        hc_pos(3:4)=hc_pos(3:4);%*1.1;
                        set(hc,'Position',hc_pos)
                        
                    end
                    
                    if cmpr1==rows_subplot && cmpr2==cols_subplot
                        %%% colorbar differences
                        subplot_handle(pos_subplot)=...
                            axes('Units','normalized','OuterPosition',[(col_act-1)*1/cols_subplot def_height_txt+(row_act-1)*(1-def_height_txt)/rows_subplot 1/cols_subplot (1-def_height_txt)/rows_subplot]);
                        set(gca,'UserData','top')
                        if no_comparemaps>1
                            map_act=squeeze(maps_grid_all_figures(no_comparemaps-1,no_comparemaps,:,:,fg));
                            h=imagesc(edgesx,edgesy,((map_act-clims_all(no_comparemaps-1,no_comparemaps,1))/(clims_all(no_comparemaps-1,no_comparemaps,2)-clims_all(no_comparemaps-1,no_comparemaps,1)))*255+(cmap_subplot(no_comparemaps-1,no_comparemaps)-1)*256+1);
                            
                        else
                            map_act=squeeze(maps_grid_all_figures(1,1,:,:,fg));
                            h=imagesc(edgesx,edgesy,((map_act-clims_all(1,1,1))/(clims_all(1,1,2)-clims_all(1,1,1)))*255+(cmap_subplot(1,1)-1)*256+1);
                            
                        end
                        set(gca,'CLim', [0 max(cmap_subplot(:))*256])
                        set(gca,'UserData','top')
                        set(h,'Visible','off')
                        set(gca,'Visible','off')
                        
                        colormap(cmap_all)
                        
                        hc=colorbar;
                        set(gca,'UserData','top')
                        set(hc,'location','southoutside');
                        labelstring=cellstr(num2str((clims_all(cmpr1-1,cmpr2,1):diff(clims_all(cmpr1-1,cmpr2,:))/4:clims_all(cmpr1-1,cmpr2,2))', '%.2g'));
                        set(hc,'XLim', [(cmap_subplot(cmpr1-1,cmpr2)-1)*256+1 (cmap_subplot(cmpr1-1,cmpr2))*256+1],...
                            'Xtick',(cmap_subplot(cmpr1-1,cmpr2)-1)*256+1:(256-1)/4 :(cmap_subplot(cmpr1-1,cmpr2))*256+1,...
                            'XtickLabel', labelstring,'FontSize',tickfontsize);
                        title(hc,{'Differences'},'fontsize',labelfontsize)
                        hc_pos=get(hc,'Position');
                        hc_pos(1)=hc_pos(1)+.01;
                        hc_pos(2)=hc_pos(2)-.03;
                        %             hc_pos(2)=hc_pos(2)+.08;
                        hc_pos(3:4)=hc_pos(3:4);%*1.1;
                        set(hc,'Position',hc_pos)
                    end
                    if cmpr1==rows_subplot-1 && cmpr2==cols_subplot-1
                        %%% colorbar p-Values
                        subplot_handle(pos_subplot) =...
                            axes('Units','normalized','OuterPosition',[(col_act-1)*1/cols_subplot def_height_txt+(row_act-1)*(1-def_height_txt)/rows_subplot 1/cols_subplot (1-def_height_txt)/rows_subplot]);
                        set(gca,'UserData','top')
                        if no_comparemaps>1
                            map_act=squeeze(maps_grid_all_figures(no_comparemaps,no_comparemaps-1,:,:,fg));
                            h=imagesc(edgesx,edgesy,((map_act-clims_all(no_comparemaps-1,no_comparemaps,1))/(clims_all(no_comparemaps,no_comparemaps-1,2)-clims_all(no_comparemaps,no_comparemaps-1,1)))*255+(cmap_subplot(no_comparemaps,no_comparemaps-1)-1)*256+1);
                        else
                            map_act=squeeze(maps_grid_all_figures(1,1,:,:,fg));
                            h=imagesc(edgesx,edgesy,((map_act-clims_all(1,1,1))/(clims_all(1,1,2)-clims_all(1,1,1)))*255+(cmap_subplot(1,1)-1)*256+1);
                            
                        end
                        set(gca,'CLim', [0 max(cmap_subplot(:))*256])
                        set(gca,'UserData','top')
                        set(h,'Visible','off')
                        set(gca,'Visible','off')
                        
                        colormap(cmap_all)
                        
                        hc=colorbar;
                        set(gca,'UserData','top')
                        set(hc,'location','southoutside');
                        labelstring=cellstr(num2str((clims_all(cmpr1+1,cmpr2,1):diff(clims_all(cmpr1+1,cmpr2,:))/4:clims_all(cmpr1+1,cmpr2,2))', '%.2g'));
                        set(hc,'XLim', [(cmap_subplot(cmpr1+1,cmpr2)-1)*256+1 (cmap_subplot(cmpr1+1,cmpr2))*256+1],...
                            'Xtick',(cmap_subplot(cmpr1+1,cmpr2)-1)*256+1:(256-1)/4 :(cmap_subplot(cmpr1+1,cmpr2))*256+1,...
                            'XtickLabel', labelstring,'FontSize',tickfontsize);
                        title(hc,{'p-Values'},'fontsize',labelfontsize)
                        hc_pos=get(hc,'Position');
                        hc_pos(1)=hc_pos(1)+.01;
                        %             hc_pos(2)=hc_pos(2)+.0;
                        hc_pos(3:4)=hc_pos(3:4);%*1.1;
                        set(hc,'Position',hc_pos)
                    end
                    
                    % Annotations and legend
                    if row_act == 1 && col_act>1
                        %                     if max([1,rows_subplot,cols_subplot])>1
                        %                         annotstr=[annotstr 'Fig. (' num2str(row_act) ',' num2str(col_act) ') '];
                        %                     end
                        annotstr=[output_plot.data_string{fg, col_act-1} '; ' annotstr ];
                        
                        
                        
                    end
                    set(subplot_handle(pos_subplot),'UserData','top')
                end
            end
            if significance_between_groups && ~isempty(output_plot.data_sign)
                
                annotstr=[annotstr 'Statistical test: ' statistical_test_type ' at significance level ' num2str(statistical_test_significance)];
                
            end
        else % If there is only one map to display
            subplot_handle=...
                axes('Units','normalized','OuterPosition',[0 def_height_txt 1 (1-def_height_txt)]);
            
            map_act=squeeze(output_plot.data{fg,1});
            if no_cols==1;
                map_act=map_act';
                
            end
            
            %%
            maxv=max(map_act(:));
            minv=min(map_act(:));
            
            range_ratio=(maxv-minv)/minv;
            
                if min(map_act(:))>=0 && range_ratio>.7
                    
                    %Normalization:
                    map_act = map_act./sum(map_act(:))./(edgesx(2)-edgesx(1))./(edgesy(2)-edgesy(1));
                    
                    clim=clim_meanmaps_norm;
                    clims = [0 clim];
                    colormap(jet)
                elseif min(map_act(:))>=0 && range_ratio<.7 && ~( clim_meanmaps_bottom<.5 && clim_meanmaps>.5)
                    
                    clims =[clim_meanmaps_bottom clim_meanmaps];
                    colormap(colormap4InteractionMaps)
                elseif min(map_act(:))>=0 && range_ratio<.7 && ( clim_meanmaps_bottom<.5 && clim_meanmaps>.5)
                    clims=[.5-max(abs(.5-clim_meanmaps_bottom),abs(.5-clim_meanmaps)) .5+max(abs(.5-clim_meanmaps_bottom),abs(.5-clim_meanmaps))];
                    colormap(colormap4InteractionMaps)
                elseif min(map_act(:))<0
                    
                    clims=[-clim clim];
                    colormap(colormap4InteractionMaps)
                else
                    clims=[0 1];
                end
               
            if ~isempty(plot_clims)
                clims=plot_clims;
                
            end
            map_act(map_act> clims(2))=clims(2);
                map_act(map_act< clims(1))=clims(1);
            
            
            imagesc(edgesx,edgesy,map_act,clims);
            
            
            no_xticks=length(get(gca,'XTick')); no_yticks=length(get(gca,'YTick'));
            
            hact=gca;
            set(gca,'YDir','normal','YtickLabelMode','Manual','YtickMode','Manual')
            axis square tight
            
            no_xticks_manual=length(ticsx);
            no_yticks_manual=length(ticsy);
            
            if no_xticks_manual==no_yticks_manual
                no_xticks=max(no_xticks,no_yticks);
                no_yticks=no_xticks;
            end
            
            tick_ratio=max(1,floor(no_xticks_manual/no_xticks));
            zeropos=find(ticsx==0);
            if ~isempty(zeropos)
                ticsx_1sthalf=ticsx(zeropos:-tick_ratio:1); ticsx_1sthalf=ticsx_1sthalf(end:-1:1);
                ticsx_2ndhalf= ticsx(zeropos+tick_ratio:tick_ratio:end);
                ticsx_new=[ticsx_1sthalf ticsx_2ndhalf];
            else
                ticsx_new=ticsx(1:tick_ratio:end);
            end
            
            tick_ratio=max(floor(no_yticks_manual/no_yticks),1);
            zeropos=find(ticsy==0);
            if ~isempty(zeropos)
                ticsy_1sthalf=ticsy(zeropos:-tick_ratio:1); ticsy_1sthalf=ticsy_1sthalf(end:-1:1);
                ticsy_2ndhalf= ticsy(zeropos+tick_ratio:tick_ratio:end);
                ticsy_new=[ticsy_1sthalf ticsy_2ndhalf];
            else
                ticsy_new=ticsy(1:tick_ratio:end);
            end
            
            ylabel(gca,options.plot_axislabelstring{2},'fontsize',labelfontsize)
            set(gca,'Xtick',ticsx_new,'XtickLabel',ticsx_new,'Ytick',ticsy_new,'YtickLabel',ticsy_new,'FontSize',tickfontsize);
            xlabel(options.plot_axislabelstring{1},'fontsize',labelfontsize)
            
            if ~isempty(plot_mode.legend{1})
                title([options.plot_titlestring{1} ': ' plot_mode.legend{1}  ])
            else
                title(options.plot_titlestring{1});
            end
            hc=colorbar;
            set(hc,'YLim', clims,'FontSize',tickfontsize,'UserData','top');
            
            
            %         B = axes('Position',get(hact,'Position'),'YLim',get(hact,'YLim'),'XLim',get(hact,'XLim'));
            %         set(B,'Xtick',ticsx_new,'Ytick',ticsy_new,'YtickLabel',[],'XtickLabel',[],'FontSize',labelfontsize)
            
            %         axis tight square
            %         uistack(B,'bottom')
            
            % Annotations and legend
            annotstr=[annotstr output_plot.data_string{fg, 1}(3:end) '; '];
            set(subplot_handle,'UserData','top')
            %         set(B,'UserData','top')
        end
        
        set(findobj(fh_raw,'Type','axes'),'Userdata','top')
        
        handle_alllpots=findobj(fh_raw,'Type','axes','-not','Tag','Colorbar');
        allpos=get(handle_alllpots,'Position');
        if ~iscell(allpos)
            allpos={allpos};
        end
        minwidth=min(cellfun(@(x) x(3),allpos));
        minheight=min(cellfun(@(x) x(4),allpos));
        for kk=1:size(handle_alllpots,1)
            
            set(handle_alllpots(kk),'Position',[allpos{kk}(1:2) minwidth minheight])
        end
        
        ax=axes('Units','Normal','Position',[-.33 .081 .85 .85],'Visible','off');
        set(get(ax,'Title'),'Visible','on','horizontalAlignment', 'center')
        %     title(['(:,:)'],'fontsize',gr_tr_fontsize);
        ax=axes('Units','Normal','Position',[.095 .075 .85 .85],'Visible','off');
        set(get(ax,'Title'),'Visible','on','horizontalAlignment', 'center')
        
        
        %% Put color background
        
        set(fh_raw,'ResizeFcn',@(one,two) resize_annotations(one,two,annotstr))
        if caption
            ha=axes('Units','normalized','OuterPosition',[0 0 1 def_height_txt],'XTick',[],'YTick',[],'xticklabel',[],'yticklabel',[]);%,'Visible','off');
            set(ha,'UserData','bottom')
            box on
            
            set(fh_raw,'UserData',annotstr);
            
            txt=text(0.01,0.02, annotstr  ,'Units','normalized','VerticalAlignment','bottom','HorizontalAlignment','left');
        end
        % The following is a hack to call 'resize_annotations'
        fgpos=get(fh_raw,'Position');
        set(fh_raw,'Position',[fgpos(1) fgpos(2) fgpos(3)*.99 fgpos(4)])
        
        %% Save a copy. As jpeg and fig.
        idSocial_saveFigure([act_method '_' fname],project_path)
        
    end
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
    
    cbh=findobj(one,'Type','axes','-and','UserData','top','-and','Tag','Colorbar');
    if size(ax1,1)==2 && ~isempty(cbh)
        ph=setxor(ax1,cbh);
        
        act_pos=get(ph,'OuterPosition');
        posy=act_pos(2);
        posx=act_pos(1);
        width=act_pos(3);
        new_height=act_pos(4)/(1-def_height_txt)*(1-def_height_txt*extnty);
        
        %         (posy-def_height_txt)/def_height_txt == (new_posy-(def_height_txt*extnty))/def_height_txt*extnty
        new_posy =  def_height_txt*extnty+(posy-def_height_txt)/(1-def_height_txt) * (1-def_height_txt*extnty);
        %         new_posy=def_height_txt*extnty+(posy-def_height_txt)*(1-def_height_txt*extnty);%/(1-def_height_txt)*extnty;
        set(ph,'OuterPosition',[posx new_posy width new_height])
        
        
        % Now the colorbar
        act_pos=get(cbh,'OuterPosition');
        act_pos_ph=get(ph,'Position');
        posy=act_pos(2);
        posx=(act_pos_ph(3))+.05;
        width=act_pos(3);
        new_height=act_pos(4)/(1-def_height_txt)*(1-def_height_txt*extnty);
        
        %         (posy-def_height_txt)/def_height_txt == (new_posy-(def_height_txt*extnty))/def_height_txt*extnty
        new_posy =  def_height_txt*extnty+(posy-def_height_txt)/(1-def_height_txt) * (1-def_height_txt*extnty);
        %         new_posy=def_height_txt*extnty+(posy-def_height_txt)*(1-def_height_txt*extnty);%/(1-def_height_txt)*extnty;
        
        set(cbh,'OuterPosition',[posx new_posy width new_height])
    else
        for k=1:size(ax1,1)
            if strcmp('Colorbar',get(ax1(1),'Tag'))
                act_pos=get(ax1(k),'OuterPosition');
                posy=act_pos(2);
                posx=act_pos(1);
                width=act_pos(3);
                new_height=act_pos(4)/(1-def_height_txt)*(1-def_height_txt*extnty);
                
                %         (posy-def_height_txt)/def_height_txt == (new_posy-(def_height_txt*extnty))/def_height_txt*extnty
                new_posy =  def_height_txt*extnty+(posy-def_height_txt)/(1-def_height_txt) * (1-def_height_txt*extnty);
                %         new_posy=def_height_txt*extnty+(posy-def_height_txt)*(1-def_height_txt*extnty);%/(1-def_height_txt)*extnty;
                
                set(ax1(k),'OuterPosition',[posx new_posy width new_height])
            end
        end
    end
    %     handle_alllpots=findobj(one,'Type','axes','-not','Tag','Colorbar');
    %     allpos=get(handle_alllpots,'Position');
    %     minwidth=min(cellfun(@(x) x(3),allpos));
    %     minheight=min(cellfun(@(x) x(4),allpos));
    %     for kk=1:size(handle_alllpots,1)
    %
    %         set(handle_alllpots(kk),'Position',[allpos{kk}(1:2) minwidth minheight])
    %     end
    
end
end
