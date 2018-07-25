function [input_data, figure_windows] = idSocial_plotUI(input_data,funcstring,act_set,control,act_fig)

if nargin < 3 || isempty(act_set)
    act_set = []; % Change this to last set!
end
if nargin < 4 || isempty(control)
    control = false;
end
act_set = act_set+0.5*control;
stat_menuH = [];
% Is act_fig figure handle?
figure_windows = [];
new_fig = true;
distribution_flag = false;

if ~isempty(act_set)
    act_data = input_data.(funcstring).results(:,floor(act_set)+1);
%     act_data = input_data.(funcstring).(['Set' num2str(floor(act_set))]).output_plot.data;
%     if iscell(act_data)
%         act_data  = cellfun(@(x) squeeze(x)',act_data,'UniformOutput',false);
%         act_data = vertcat(act_data{:});
%     end
    pm = input_data.(funcstring).(['Set' num2str(floor(act_set))]).plot_mode;
    dm = pm.display_mode;
    % For distributions:
    distr_statisics = pm.statistics{1};
    distribution_flag = ~isempty(strfind(lower(distr_statisics),'hist'));
    distr_statistics = strrep(lower(distr_statisics), 'hist+','');

end


sign_handle = [];
usdat = [];
% act_fig is there and is a handle: Plot onto existing figure, get number
% from 'UserData'
if nargin > 4 && ~isempty(act_fig) && IsGraphicHandle(act_fig)
    if isprop(act_fig,'UserData')
        usdat = get(act_fig,'UserData');
    else
        usdat = [];
    end
    if ~isempty(usdat) && isstruct(usdat) && isfield(usdat,'act_fig')
        figure_windows = act_fig;
        act_fig = usdat.act_fig;
        figure(figure_windows)
        figure_windows=removeCallbacksToFigure(figure_windows);
        [figure_windows,stat_menuH] = addCallbacksToFigure(figure_windows);
    end
end

% act_fig is a number, but figure_windows is still empty: Figure should
% exist, but we have to find the handle either in input_data or
% "everywhere", or we have to load it.
if nargin > 4 && isnumeric(act_fig) && isempty(figure_windows)
    % Try find figure within the open plot result figures:
    allFigs = get(0,'Children');
    plotFigs = findobj(allFigs,'-regexp','Name',['Plot Results ' num2str(act_fig)]);
    if any(arrayfun(@(x) ~isempty(x),plotFigs) & IsGraphicHandle(plotFigs))
        figure_windows = plotFigs;
        thisFig = false(1,numel(plotFigs));
        for k3 = 1:numel(plotFigs)
            tud = get(plotFigs(k3),'UserData');
            if isfield(tud,'date') && ~isempty(tud.date) && ischar(tud.date) && ...
                    strcmpi(tud.date,input_data.(funcstring).(['Fig' num2str(act_fig)]).date)
                thisFig(k3)=true;
            end
        end
        if sum(thisFig)==1
            figure_windows = plotFigs(thisFig);
            figure_windows=removeCallbacksToFigure(figure_windows);
            [figure_windows,stat_menuH] = addCallbacksToFigure(figure_windows);
            % Update handle to be sure:
%             input_data.(funcstring).(['Fig' num2str(act_fig)]).handle = figure_windows;
        else
            figure_windows = [];
        end
        
    end
    % If the figure could not be found try loading it
    if isempty(figure_windows) && isfield(input_data.(funcstring),['Fig' num2str(act_fig)]) && isfield(input_data.(funcstring).(['Fig' num2str(act_fig)]),'path') && ...
            exist(input_data.(funcstring).(['Fig' num2str(act_fig)]).path,'file')==2
        figure_windows =  openfig(input_data.(funcstring).(['Fig' num2str(act_fig)]).path);
        ud4 = get(figure_windows,'UserData');
        dm = ud4.display_mode;
        input_data.(funcstring).(['Fig' num2str(act_fig)]).act_set = ud4.act_set;
        control = mod(ud4.act_set,1)==0.5; 
        act_sets_ud = floor(ud4.act_set);
        if isempty(act_set)
            act_set = act_sets_ud;
            act_data = input_data.(funcstring).results(:,floor(act_set)+1);
%             pm = input_data.(funcstring).(['Set' num2str(floor(act_set))]).plot_mode;
%            
%             % For distributions:
%             distr_statisics = pm.statistics{1};
%             distribution_flag = ~isempty(strfind(lower(distr_statisics),'hist'));
%             distr_statistics = strrep(lower(distr_statisics), 'hist+','');
        end
        
        ud4.save_path = input_data.(funcstring).(['Fig' num2str(act_fig)]).path;
        set(figure_windows,'UserData',ud4);
        
        distr_statisics_cell = cell(numel(act_sets_ud),1);
        for k4 = 1:numel(act_sets_ud)
            pm = input_data.(funcstring).(['Set' num2str(act_sets_ud(k4))]).plot_mode;
            % For distributions:
            distr_statisics_cell{k4} = pm.statistics{1};
        end
        if numel(distr_statisics_cell)>0 && all(cellfun(@(x) strcmp(distr_statisics_cell{1},x),distr_statisics_cell))
            distribution_flag = ~isempty(strfind(lower(distr_statisics_cell{1}),'hist'));
            distr_statistics = strrep(lower(distr_statisics_cell{1}), 'hist+','');
        end
        [figure_windows,stat_menuH] = addCallbacksToFigure(figure_windows);
        
        
%         udFig = get(gcf,'UserData');
%         act_set = udFig.act_set;
%         act_data = input_data.(funcstring).results(:,floor(act_set)+1);
% 
%         pm = input_data.(funcstring).(['Set' num2str(floor(act_set))]).plot_mode;
%         dm = pm.display_mode;
%         % For distributions:
%         distr_statisics = pm.statistics{1};
%         distribution_flag = ~isempty(strfind(lower(distr_statisics),'hist'));
%         distr_statistics = strrep(lower(distr_statisics), 'hist+','');
        
        new_fig = false;
    end
    
    if isempty(figure_windows)
        error([mfilename ': Sorry, there was something wrong. The previous figure could not be found.'])
    end
    
end

% Finally, if act_fig does not exist or is empty: Create new figure; get number of
% new figure from input_data.(funcstring)
if nargin < 5 || isempty(act_fig)
%     fn = fieldnames(input_data.(funcstring));
%     maxFig = 1;
%     for k=1:size(fn,1)
%         if ~isempty(strfind(fn{k},'Fig'))
%             maxFig = max(maxFig,str2double(strrep(fn{k},'Fig','')));
%         end
%     end
%     act_fig = maxFig+1;
    % Note: if act_set is empty, simply a new figure will be opened.
    [figure_windows, input_data, act_fig, stat_menuH] = newFigure(input_data,funcstring,act_set);
end


user_data = get(figure_windows,'UserData');
% If figure has no plot yet
if ~isfield(user_data,'display_mode') || isempty(user_data.display_mode)
    user_data.display_mode = pm.display_mode;
    set(figure_windows,'UserData',user_data)
end

[~, temp_savepath]= idSocial_recursiveGetOptionsFromOptionsCell(input_data.options,'temp_savepath');

input_data.(funcstring).(['Fig' num2str(act_fig)]).path = [temp_savepath funcstring '_figure' num2str(act_fig) '.fig'];



feature('DefaultCharacterSet','UTF-8');



[~,project_name]=idSocial_recursiveGetOptionsFromOptionsCell(input_data.options,'project_name');

if ~isempty(act_set)
    if ~control
        dat = act_data{3};
        plotResultsIndivData = act_data{4};
        
    elseif control && size(act_data,1)>=10 && ~isempty(act_data{10})
        dat = act_data{10};
        plotResultsIndivData = act_data{11};
    else
        warning([mfilename ': Could not find control data'])
        return
    end
    
    labels = act_data(9);
    
    xticks = act_data{8}{1};
    yticks = act_data{8}{2};
    
    xticklabel = act_data{9}{1};
    yticklabel = act_data{9}{2};
    
    if iscell(xticklabel) && numel(xticklabel)==1; xticklabel = xticklabel{1}; end
    if iscell(xticks); xticks = xticks{1}; end
    
    if iscell(yticklabel)&& numel(yticklabel)==1; yticklabel = yticklabel{1}; end
    if iscell(yticks); yticks = yticks{1}; end
    
    
    xaxislabel = act_data{7}{1};
    yaxislabel = act_data{7}{2};
end









%% For statistics
if isempty(figure_windows)
    no_plots = 0;
    lineH = NaN(1,1);
else
    [lineH, no_plots] = update_LineHandle(figure_windows(end));

end


%%% Check if new data is compatible to already existing plot
if ~isempty(act_set) && ~isempty(lineH) && numel(lineH)>1 &&  ishghandle(lineH(end-1))
    ax = get(lineH(end-1),'Parent');
    prev_xlabel = get(ax,'XLabel');
    prev_ylabel = get(ax,'YLabel');
    user_data = get(figure_windows,'UserData');
    if ~strcmpi(xaxislabel,get(prev_xlabel,'String')) || ~strcmpi(yaxislabel,get(prev_ylabel,'String')) || ~strcmpi(user_data.display_mode,pm.display_mode)
        warning([mfilename ': Sorry, new data does not fit the data already plotted.'])
        return;
    end
end
%%

errH = [];

% For maps:
comp_set_idx = -1;
Xctrs = [];
Yctrs = [];

figure_windows = figure_windows(ishandle(figure_windows));



% cmap = get(figure_windows(end),'Colormap');
colororder = jet(no_plots+1);

control_flag = false(no_plots+1,1);
linewidth = NaN(no_plots+1,1);
linewidth_DistrStats = NaN(no_plots+1,1);
linecolor = NaN(no_plots+1,3);
xdataCell = cell(no_plots+1,1);
ydataCell = cell(no_plots+1,1);
errH = NaN(no_plots+1,1);
distrStatsH = NaN(no_plots+1,1);
errorbar_type = repmat({'None'},[no_plots+1,1]);

for k = 1:no_plots
    if verLessThan('matlab','8.4')
        set(lineH(k),'Color',colororder(k,:));
    end
    linewidth(k) = get(lineH(k),'LineWidth');
    linewidth_DistrStats(k) = linewidth(k);
    linecolor(k,:) = get(lineH(k),'Color');
    set(lineH(k),'DeleteFcn',{@LineDelete,k});
    xdataCell{k} = get(lineH(k),'XData');
    ydataCell{k} = get(lineH(k),'YData');
end

marked_lines = zeros(1,no_plots+1);
marked_DistrStats = zeros(1,no_plots+1);

marked_points = zeros(max(cellfun(@(x) numel(x),xdataCell)),no_plots+1);
markh = NaN(max(cellfun(@(x) numel(x),xdataCell)),no_plots+1);

no_possiblemarks = inf;
 ylim = get(gca,'YLim');
%%
if ~isempty(act_set) && (isempty(usdat) || ~isempty(usdat) && isstruct(usdat) && isfield(usdat,'act_set') && ~any((usdat.act_set)==act_set)) && new_fig
%     allfigs = findall(0,'-regexp','Name','Plot Results');
%     for af = 1:numel(allfigs)
%         allfigs(af)=removeCallbacksToFigure(allfigs(af));
%     end
    if  strcmpi(dm,'dynamicMapPolar')
 
        if strcmpi(pm.statistics,'positive_ratio') || ...
                (strcmpi(pm.statistics,'pool') && strcmpi(pm.data,'positive_ratio'))
            limval = max(1-min(dat(:)), max(dat(:)));
            clims=[1-limval limval];
        else
            limval = max(abs(dat(:)));
            clims=[-limval  limval];
        end
        
        [h,ax_lines,Xctrs,Yctrs]= ...
            idSocial_imagescPolar(xticks,yticks,dat',clims,[12 12]);
        info.project_name = project_name;
        info.funcstring = funcstring;
        info.act_set = act_set;%+control*0.5;
        info.control = control;
        set(h,'UserData',info)
        
        set(gca,'YDir','normal')
        set(h, 'EdgeColor', 'none');
        axis equal  tight;% off
        hold on
        xlabel(xaxislabel)
        ylabel(yaxislabel)
        set(gca,'XTick',[-xticks(end:-1:2) xticks],'XTickLabel',[-xticks(end:-1:2) xticks])
        set(gca,'YTick',[-xticks(end:-1:2) xticks],'YTickLabel',[-xticks(end:-1:2) xticks])
        set(gca, 'XTickMode', 'auto', 'XTickLabelMode', 'auto')
        set(gca, 'YTickMode', 'auto', 'YTickLabelMode', 'auto')
        colormap4InteractionMapsRedBlue = [];
        load('colormap4InteractionMapsRedBlue.mat')
        colormap(colormap4InteractionMapsRedBlue)
        colorbar

        idSocial_improveGraphics(figure_windows(end))
        
        
        %         input_data.(funcstring).(['Fig' num2str(act_fig)]).handle = figure_windows(end);
        if verLessThan('matlab','8.2')
            
            saveas(figure_windows(end), input_data.(funcstring).(['Fig' num2str(act_fig)]).path,'fig');
        else
            savefig(figure_windows(end), input_data.(funcstring).(['Fig' num2str(act_fig)]).path);
        end
        
    elseif  strcmpi(dm,'MapPolar')
        

        limval = max(abs(dat(:)));
        clims=[0  limval];
        
        [h,ax_lines,Xctrs,Yctrs]= ...
            idSocial_imagescPolar(xticks,yticks,dat',clims,[12 12]);
        info.project_name = project_name;
        info.funcstring = funcstring;
        info.act_set = act_set;%+control*0.5;
        info.control = control;
        set(h,'UserData',info)
        
        set(gca,'YDir','normal')
        set(h, 'EdgeColor', 'none');
        axis equal  tight;% off
        hold on
        xlabel(xaxislabel)
        ylabel(yaxislabel)
        set(gca,'XTick',[-xticks(end:-1:2) xticks],'XTickLabel',[-xticks(end:-1:2) xticks])
        set(gca,'YTick',[-xticks(end:-1:2) xticks],'YTickLabel',[-xticks(end:-1:2) xticks])
        set(gca, 'XTickMode', 'auto', 'XTickLabelMode', 'auto')
        set(gca, 'YTickMode', 'auto', 'YTickLabelMode', 'auto')
        colormap(jet)
        colorbar

        idSocial_improveGraphics(figure_windows(end))
        
        
        
%         input_data.(funcstring).(['Fig' num2str(act_fig)]).handle = figure_windows(end);
        if verLessThan('matlab','8.2')
            saveas(figure_windows(end), input_data.(funcstring).(['Fig' num2str(act_fig)]).path,'fig');
        else
            savefig(figure_windows(end), input_data.(funcstring).(['Fig' num2str(act_fig)]).path);
        end
    else
  
        figure(figure_windows(end))
     
        
        hold on
        try
            ph=plot(xticks,dat,'-o');
            info.project_name = project_name;
            info.funcstring = funcstring;
            info.act_set = act_set;%+control*0.5;
            info.control = control;
            
            set(ph,'UserData',info)
            
            if verLessThan('matlab','8.4')
                set(ph,'Color',colororder(no_plots+1,:));
            end
            linewidth(end) =  get(ph,'LineWidth');
            linewidth_DistrStats(end) =  get(ph,'LineWidth');
            linecolor(end,:) = get(ph,'Color');
            set(ph,'DeleteFcn',{@LineDelete,no_plots+1}); %%
%             UNCOMMENT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            control_flag(end) = control & size(act_data,1)>=10;
            
            xdataCell{end} = get(ph,'XData');
            ydataCell{end} = get(ph,'YData');
            
            
            
            %         line_ud = get(ph,'UserData');
            %         line_ud.listenerH1 = addlistener(ph, 'LineWidth', 'PostSet', @(hAxes, eventData) setLineWidth(hAxes, eventData,no_plots+1,lineH(no_plots+1)));
            %         line_ud.listenerH2 = addlistener(ph, 'Color', 'PostSet', @(hAxes, eventData) setColor(hAxes, eventData,no_plots+1,lineH(no_plots+1)));
            %         set(ph,'UserData',line_ud);
            lineH(end) = ph;
            
            errorbar_type{end} = 'None';
            
            
            xlabel(xaxislabel)
            ylabel(yaxislabel)
            set(gca,'XTick',xticks,'XTickLabel',xticklabel)
            idSocial_improveGraphics(figure_windows(end))
            ylim = get(gca,'YLim');
        catch
            close(figure_windows(end))
            warndlg(['Sorry, I cannot plot the data.'])
            
        end
        
    
        
        %         input_data.(funcstring).(['Fig' num2str(act_fig)]).handle = figure_windows(end);
        
        %         ud3 = get(input_data.(funcstring).(['Fig' num2str(act_fig)]).handle,'UserData');
        %         ud3.save_path =  [temp_savepath funcstring '_figure' num2str(act_fig) '.fig'];
        %         set(input_data.(funcstring).(['Fig' num2str(act_fig)]).handle,'UserData',ud3)
        ud3 = get(figure_windows(end),'UserData');
        ud3.save_path =  [temp_savepath funcstring '_figure' num2str(act_fig) '.fig'];
        set(figure_windows(end),'UserData',ud3)
        
    end
    
    udat = get(figure_windows(end),'UserData');
    if ~isstruct(udat); udat = struct; end
    if isfield(udat,'act_set') && ~isempty(udat.act_set)
        udat.act_set = unique([udat.act_set act_set]); %% Add 0.5 if control to identify control data later
        input_data.(funcstring).(['Fig' num2str(act_fig)]).act_set = udat.act_set;
    else
        udat.act_set =  act_set;
        input_data.(funcstring).(['Fig' num2str(act_fig)]).act_set = udat.act_set;
    end
%     if isfield(udat,'date') && ~isempty(udat.date)
%         input_data.(funcstring).(['Fig' num2str(act_fig)]).date = udat.date;
%     end
    udat.act_fig = act_fig;
    udat.act_func = funcstring;
    set(figure_windows(end),'UserData',udat);
%     allfigs = findall(0,'-regexp','Name','Plot Results');
%     for af = 1:numel(allfigs)
%         allfigs(af)=addCallbacksToFigure(allfigs(af));
%     end
end


    function setStatTest(src,event,act_options,fh)
        
        out_opts = idSocialUI_propertyTable([],act_options);
        if ~isempty(out_opts)
            [lineH,no_plots] = update_LineHandle(fh);
            out_opts = handleComboMenus(out_opts);
            if ~isempty(out_opts)
                
                if strcmpi(dm,'dynamicMapPolar') || strcmpi(dm,'MapPolar')
                else
                    act_info = get(lineH(~isnan(lineH)),'UserData');
                    if ~iscell(act_info)
                        adum = act_info;
                        act_info = cell(numel(adum));
                        for ac = 1:numel(adum)
                            act_info{ac} = adum(ac);
                        end
                    end
                    
                    if any(marked_points(:))
                        [~,point_idces] = sort(marked_points(:),'descend');
                        point_idces = point_idces(1:sum(marked_points(:)>0));
                        lidx = ceil(point_idces/size(marked_points,1));
                        pidx = mod(point_idces-1,size(marked_points,1))+1;
                    elseif any(marked_DistrStats(:))
                        [~,point_idces] = sort(marked_DistrStats(:),'descend');
                        point_idces = point_idces(1:sum(marked_DistrStats(:)>0));
                        lidx = ceil(point_idces/size(marked_DistrStats,1));
                        
                    elseif any(marked_lines)
                        [~, line_idces] = sort(marked_lines,'descend');
                    else
                        line_idces = 1:size(act_info,1);
                        
                    end
                end
                
                if strcmpi(dm,'dynamicMapPolar') || strcmpi(dm,'MapPolar')
                    no_sets = size(input_data.(funcstring).results,2)-1;
                    display_mode = dm;
                    val1 = input_data.(funcstring).results{4,floor(act_set)+1};
                    if comp_set_idx>0 && comp_set_idx<no_sets+1
                        val2 = input_data.(funcstring).results{4,comp_set_idx+1};
                    elseif comp_set_idx==no_sets+1
                        if  size(act_data,1)>=10 && ~isempty(act_data{10})
                            val2 = input_data.(funcstring).results{11,floor(act_set)+1};
                        else
                            warning([mfilename ': Could not find control data'])
                            return
                        end
                    end
                    val = NaN(2,max(size(val1,1),size(val1,2)),size(val1,2),size(val1,3));
                    val(1,1:size(val1,1),:,:) =val1;
                    val(2,1:size(val2,1),:,:) =val2;
                    lineHReorder = [];
                    
                else
                    display_mode = 'plot2d';
                    if any(marked_points(:)) || any(marked_DistrStats(:))
                        lineHReorder = lineH;
                        datasetIndiv = cell(size(lidx,1),1);
                        xticks_act = cell(size(lidx,1),1);
                        count = 1;
                        for lh = lidx'
                            
                            if any(marked_points(:))
                                if ~control_flag(lidx(count))
                                    datasetIndiv{count} = input_data.(act_info{lh}.funcstring).results{4,floor(act_info{lh}.act_set)+1}(:,pidx(count));
                                else
                                    datasetIndiv{count} = input_data.(act_info{lh}.funcstring).results{11,floor(act_info{lh}.act_set)+1}(:,pidx(count));
                                    
                                end
                                xtdum = get(lineH(lh),'XData');
                                xticks_act{count} = xtdum(pidx(count));
                            elseif any(marked_DistrStats(:))
                                datasetIndiv{count} = input_data.(act_info{lh}.funcstring).results{6,floor(act_info{lh}.act_set)+1};
                                xtdum = get(distrStatsH(lh),'XData');
                                xticks_act{count} = xtdum(1);
                            end
                            
                            
                            count = count + 1;
                        end
                        no_max_dpts = max(cellfun(@(x) size(x,1),datasetIndiv));
                        val = NaN(size(lidx,1),no_max_dpts,1);
                        for lh = 1:size(lidx,1)
                            
                            val(lh,1:size(datasetIndiv{lh},1),1) = datasetIndiv{lh};
                        end
                        display_mode = [xticks_act{:}];
                        
                    else% any(marked_lines)
                        no_markedlines = sum(marked_lines>0);
                        if no_markedlines == 0
                            no_markedlines = size(act_info,1);
                        end
                        datasetIndiv = cell(no_markedlines,1);
                        xticks_act = cell(no_markedlines,1);
                        lineHReorder = lineH(line_idces); %Re-order for idSocial_auxiliaries_plotSignificance
                        line_idces = line_idces(1:no_markedlines);
                        lineHReorder = lineHReorder(1:no_markedlines);
                        for lh = 1:no_markedlines;% line_idces
                            if ~control_flag(lh)
                                datasetIndiv{lh} = input_data.(act_info{line_idces(lh)}.funcstring).results{4,floor(act_info{line_idces(lh)}.act_set)+1};
                            else
                                datasetIndiv{lh} = input_data.(act_info{line_idces(lh)}.funcstring).results{11,floor(act_info{line_idces(lh)}.act_set)+1};
                                
                            end
                            xticks_act{lh} = get(lineH(line_idces(lh)),'XData');
                            
                        end
                        min_x = min(cellfun(@(x) min(x),xticks_act));
                        max_x = max(cellfun(@(x) max(x),xticks_act));
                        %                     if numel(xticks_act)>1
                        sx = cellfun(@(x) numel(x),xticks_act);
                        if all(sx>1)
                            dx = cellfun(@(x) x(2)-x(1),xticks_act);
                            
                            if ~all(dx==dx(1))
                                error([mfilename ': Data does not coincide in xtick spacing.'])
                            end
                            no_parts = length( min_x:dx(1):max_x);
                        else
                            no_parts = 1;
                        end
                        no_max_dpts = max(cellfun(@(x) size(x,1),datasetIndiv));
                        
                        val = NaN(numel(line_idces),no_max_dpts,no_parts);
                        count = 1;
                        for lh = 1:no_markedlines;% line_idces
                            start_idx = find(xticks_act{lh}(1)==(min_x:max_x));
                            val(count,1:size(datasetIndiv{lh},1),start_idx:start_idx+no_parts-1) = datasetIndiv{lh};
                            count =  count + 1;
                        end
                        for al = 1:numel(lineH)
                            if ishghandle(lineH(al))
                                ua = get(lineH(al),'UserData');
                                ua.LineSelected = 0;
                                set(lineH(al),'UserData',ua);
                            end
                        end
                    end
                end
                
                if ~isempty(out_opts)
                    if ~isempty(val) && ~any(size(val)==1)
                        if any(any(marked_lines)) || any(any(marked_points)) || any(any(marked_DistrStats))
                            out_opts.one_way_only = true;
                        else
                            out_opts.one_way_only = false;
                        end
                        
                        [H, p, CI, stat, sampstat]= ...
                            idSocial_auxiliaries_statisticalSignificance(val,get(src,'Label'),out_opts);
                        if iscell(out_opts.tail(end)) && iscell(out_opts.tail{end})
                            tail =  out_opts.tail{end}{1};
                        else
                            tail =  out_opts.tail{1};
                        end
                        
                        
                        sign_handle_temp = idSocial_auxiliaries_plotSignificance(lineHReorder(~isnan(lineHReorder)),p,out_opts.alpha,tail,display_mode,get(src,'Label'),Xctrs,Yctrs);
                        sign_handle = [sign_handle sign_handle_temp];
                        
                        %                     delete(stat_menuH{14})
                        %                     warning([mfilename ': The error lies here...line 589!'])
                        set(stat_menuH{14},'Callback',{@remSignMarkers,sign_handle});
                        if ~strcmpi(dm,'dynamicMapPolar') && ~strcmpi(dm,'MapPolar')
                            for ml=1:numel(lineH)
                                if ~isnan(lineH(ml)) && ishghandle(lineH(ml))
                                    set(lineH(ml), 'LineWidth', linewidth(ml));
                                end
                            end
                        end
                        clearSelection
                    end
                end
            end
        end
    end

    function [lineH2,no_plots_out] = update_LineHandle(fh)
        lineHT = findall(fh,'Type','Line');
        check_idx = false(1,size(lineHT,1));
        for k2=1:size(lineHT,1)
            udt = get(lineHT(k2),'UserData');
            if isempty(udt) || ~isstruct(udt) || ~isfield(udt,'project_name')
                if ~isempty(lineHT)
                    check_idx(k2)=true;
                end
            end
        end
        lineHT(check_idx)=[];
        no_plots_out = numel(lineHT);
        if ~isempty(lineHT)
            lineH2 = NaN(1,size(lineHT,1)+1);
            lineH2(1,1:size(lineHT,1)) = lineHT;
        else
            lineH2 = NaN(1,1);
        end
    end

    function setLineSelect(src,event)
        plotedit off
        datacursormode off
       
        clearSelection
        
        for k2 = 1:no_plots+1
            if ~isnan(lineH(k2)) && ishghandle(lineH(k2))
                set(lineH(k2),'ButtonDownFcn',{@LineSelected});
            end
            
        end
        
    end


    function setDistrStatsSelect(src,event)
        plotedit off
        datacursormode off
        
        clearSelection
        
        for k2 = 1:no_plots+1
            if ~isnan(distrStatsH(k2)) && ishghandle(distrStatsH(k2))
             set(distrStatsH(k2),'ButtonDownFcn',{@DistrStatsSelected,k2});
            end
        end
        
    end

    function setPointSelect(~,event)
        plotedit off
        datacursormode off
      
        
        clearSelection
        
        for k2 = 1:no_plots+1
            if ~isnan(lineH(k2)) && ishghandle(lineH(k2))
                set(lineH(k2),'ButtonDownFcn',{@PointSelected,k2});
            end
        end
        
    end

    function clearSelection(src,event)
        
        for mp = 1:size(marked_points,1)
            for mp2 = 1:size(marked_points,2)
                if ~isempty(markh) && ishandle(markh(mp,mp2))
                    delete(markh(mp,mp2));
                    marked_points(mp,mp2) = 0;
                end
            end
        end
        
        marked_lines = zeros(size(marked_lines));
        
        for ml=find(~marked_lines)
            if exist('lineH','var')==1 && ishghandle(lineH(ml)) && ishghandle(lineH(ml)) % && isgraphics(lineH(ml))
                if ~isnan(lineH(ml)) && ishghandle(lineH(ml))
                    set(lineH(ml), 'LineWidth', linewidth(ml));
                end
            end
        end
        marked_DistrStats = zeros(size(marked_DistrStats));
        for ml=find(~marked_DistrStats)
            if exist('distrStatsH','var')==1 && ishandle(distrStatsH(ml)) && ishghandle(distrStatsH(ml)) %%&& isgraphics(distrStatsH(ml))
                set(distrStatsH(ml), 'LineWidth', linewidth_DistrStats(ml));
            end
        end
        
    end

    function clearSetSelection(src,event)
        no_sets = size(input_data.(funcstring).results,2)-1;

        for as1 = 1:size(stat_menuH{7},2)
            if ~isnan(stat_menuH{7}(as1)) && as1<=no_sets
                set(stat_menuH{7}(as1),'Label',['Set ' num2str(as1)]);
            elseif ~isnan(stat_menuH{7}(as1)) && as1>no_sets
                set(stat_menuH{7}(as1),'Label','Control');
            end
        end
        set(stat_menuH{8},'Enable','off')
    end

    function setSetSelect(src,event,comp_set)
        no_sets = size(input_data.(funcstring).results,2)-1;

        comp_set_idx = comp_set;
        set(stat_menuH{8},'Enable','on')
        if comp_set<=no_sets
            set(stat_menuH{7}(comp_set),'Label',['Set ' num2str(comp_set) ' ' char(hex2dec('2713'))]);
            
        else
            set(stat_menuH{7}(comp_set),'Label',['Control ' char(hex2dec('2713'))]);
        end
        
        for as1 = 1:size(stat_menuH{7},2)
            if as1~=comp_set && ~isnan(stat_menuH{7}(as1)) && as1<=no_sets
                set(stat_menuH{7}(as1),'Label',['Set ' num2str(as1)]);
            elseif as1~=comp_set && as1>no_sets && ~isnan(stat_menuH{7}(as1))
                set(stat_menuH{7}(as1),'Label','Control');
            end
        end
        
    end

    function addSTD(src,event)
        lineH = update_LineHandle(figure_windows);
        act_info = get(lineH(~isnan(lineH)),'UserData');
        if ~iscell(act_info)
            adum = act_info;
            act_info = cell(numel(adum));
            for ac = 1:numel(adum)
                act_info{ac} = adum(ac);
            end
        end
        for lh = 1:size(act_info,1)
            if ~isempty(errH(lh)) && ishandle(errH(lh))
                delete(errH(lh));
                errH(lh) = NaN;
            end
            if ~strcmpi(errorbar_type{lh},'std');
                if ~control_flag(lh)
                    dataset = input_data.(act_info{lh}.funcstring).results{3,floor(act_info{lh}.act_set)+1};
                    datasetIndiv = input_data.(act_info{lh}.funcstring).results{4,floor(act_info{lh}.act_set)+1};
                else
                    dataset = input_data.(act_info{lh}.funcstring).results{10,floor(act_info{lh}.act_set)+1};
                    datasetIndiv = input_data.(act_info{lh}.funcstring).results{11,floor(act_info{lh}.act_set)+1};
                end
                xticks_act = input_data.(act_info{lh}.funcstring).results{8,floor(act_info{lh}.act_set)+1};
                if iscell(xticks_act); xticks_act = xticks_act{1}; end
                errH(lh) = errorbar(xticks_act,dataset ,nanstd(datasetIndiv),'Color',get(lineH(lh),'Color'));
                uistack(errH(lh),'bottom')
                errorbar_type{lh} = 'std';
            end
            if isnan(errH(lh))
                errorbar_type{lh} = 'None';
            end
            
        end
    end

    function addSEM(src,event)
        lineH = update_LineHandle(figure_windows);
        act_info = get(lineH(~isnan(lineH)),'UserData');
        if ~iscell(act_info)
            adum = act_info;
            act_info = cell(numel(adum));
            for ac = 1:numel(adum)
                act_info{ac} = adum(ac);
            end
        end
        for lh = 1:size(act_info,1)
            if numel(errH)>=lh && ~isempty(errH(lh)) && ishandle(errH(lh))
                delete(errH(lh));
                errH(lh) = NaN;
                
            end
            
            if ~strcmpi(errorbar_type{lh},'sem');
                if ~act_info{lh}.control
                    dataset = input_data.(act_info{lh}.funcstring).results{3,floor(act_info{lh}.act_set)+1};
                    datasetIndiv = input_data.(act_info{lh}.funcstring).results{4,floor(act_info{lh}.act_set)+1};
                else
                    dataset = input_data.(act_info{lh}.funcstring).results{10,floor(act_info{lh}.act_set)+1};
                    datasetIndiv = input_data.(act_info{lh}.funcstring).results{11,floor(act_info{lh}.act_set)+1};
                end
                xticks_act = input_data.(act_info{lh}.funcstring).results{8,floor(act_info{lh}.act_set)+1};
                if iscell(xticks_act); xticks_act = xticks_act{1}; end
                errH(lh) = errorbar(xticks_act,dataset ,nanstd(datasetIndiv)./sqrt(sum(~isnan(datasetIndiv))),'Color',get(lineH(lh),'Color'));
                uistack(errH(lh),'bottom')
                errorbar_type{lh} = 'sem';
            end
            if isnan(errH(lh))
                errorbar_type{lh} = 'None';
            end
            
        end
    end

    function addDistrStats(src,event)
        
        lineH = update_LineHandle(figure_windows);
        act_info = get(lineH(~isnan(lineH)),'UserData');
        if ~iscell(act_info)
            adum = act_info;
            act_info = cell(numel(adum));
            for ac = 1:numel(adum)
                act_info{ac} = adum(ac);
            end
        end
        for lh = 1:size(act_info,1)
            
            dataset = input_data.(act_info{lh}.funcstring).results{5,floor(act_info{lh}.act_set)+1};
            %                 datasetIndiv = input_data.(act_info{lh}.funcstring).results{6,floor(act_info{lh}.act_set)+1};
            distrStatsH(lh) = plot([dataset dataset], ylim,'Color',get(lineH(lh),'Color'));
            set(distrStatsH(lh),'DeleteFcn',{@distrStatDelete}); %,no_plots+1
            uistack(distrStatsH(lh),'bottom')
            
        end
        set(stat_menuH{12},'Enable','on');
    end


    function remError(src,event)
        lineH = update_LineHandle(figure_windows);
        act_info = get(lineH(~isnan(lineH)),'UserData');
        if ~iscell(act_info)
            adum = act_info;
            act_info = cell(numel(adum));
            for ac = 1:numel(adum)
                act_info{ac} = adum(ac);
            end
        end
        for lh = 1:size(act_info,1)
            if ~isempty(errH(lh)) && ishandle(errH(lh))
                delete(errH(lh));
                errH(lh) = NaN;
                
            end
            if ~isempty(distrStatsH) && ishandle(distrStatsH(lh)) && ishghandle(distrStatsH(lh)) %%&& isgraphics(distrStatsH(lh))
                delete(distrStatsH(lh))
                distrStatsH(lh)=NaN;
                set(stat_menuH{12},'Enable','off');
            end
            if isnan(errH(lh))
                errorbar_type{lh} = 'None';
            end
            
        end
    end


    function addStatSign(src,event)
        errorbar(xticks,dat,nanstd(plotResultsIndivData)./sqrt(sum(~isnan(plotResultsIndivData))),'Color',get(ph,'Color'))
    end

    function LineSelected(ObjectH, EventData,act_plot)
        
        act_plot = find(eq(ObjectH,lineH));
        if ~isempty(act_plot)
            lineH = update_LineHandle(figure_windows);
            marked_lines = zeros(1,numel(lineH));
            
            for al = 1:numel(lineH)
                if ishghandle(lineH(al))
                    ua = get(lineH(al),'UserData');
                    if isfield(ua,'LineSelected') && isnumeric(ua.LineSelected)
                        marked_lines(al) = ua.LineSelected;
                    else
                        ua.LineSelected = 0;
                        set(lineH(al),'UserData',ua);
                    end
                end
            end
            if marked_lines(act_plot)>0
                %             for hc = 1:numel(H)
                %                 set(lineH(hc), 'LineWidth', linewidtlineH(hc));
                %             end
                set(lineH(act_plot), 'LineWidth', linewidth(act_plot));
                marked_lines(act_plot) = 0;
                ua = get(lineH(act_plot),'UserData');
                ua.LineSelected = 0;
                set(lineH(act_plot),'UserData',ua);
            else
                marked_lines(marked_lines>0) = marked_lines(marked_lines>0) + 1;
                marked_lines(act_plot) = 1;
                marked_lines(marked_lines>no_possiblemarks) = 0;
                
                ua = get(lineH(act_plot),'UserData');
                ua.LineSelected = 1;
                set(lineH(act_plot),'UserData',ua);
                for al = 1:numel(lineH)
                    if ishghandle(lineH(al)) && al~=act_plot
                        ua = get(lineH(al),'UserData');
                        if isfield(ua,'LineSelected') && isnumeric(ua.LineSelected) && ua.LineSelected>0
                            ua.LineSelected = ua.LineSelected+1;
                            if ua.LineSelected>no_possiblemarks
                                no_possiblemarks = 0;
                            end
                            set(lineH(al),'UserData',ua);
                        end
                    end
                end
                
                for ml=find(~marked_lines)
                    if ~isnan(lineH(ml)) && ishghandle(lineH(ml))
                        set(lineH(ml), 'LineWidth', linewidth(ml));
                    end
                end
                for ml=find(marked_lines)
                    if ~isnan(lineH(ml)) && ishghandle(lineH(ml))
                        set(lineH(ml), 'LineWidth', linewidth(ml)+2);
                    end
                end
            end
        end
    end

    function DistrStatsSelected(ObjectH, EventData,act_plot)
        
        
        
        if marked_DistrStats(act_plot)>0
            %             for hc = 1:numel(H)
            %                 set(lineH(hc), 'LineWidth', linewidtlineH(hc));
            %             end
            set(distrStatsH(act_plot), 'LineWidth', linewidth_DistrStats(act_plot));
            marked_DistrStats(act_plot) = 0;
            
        else
            marked_DistrStats(marked_DistrStats>0) = marked_DistrStats(marked_DistrStats>0) + 1;
            marked_DistrStats(act_plot) = 1;
            marked_DistrStats(marked_DistrStats>no_possiblemarks) = 0;
            for ml=find(~marked_DistrStats&~isnan(distrStatsH'))
                set(distrStatsH(ml), 'LineWidth', linewidth_DistrStats(ml));
            end
            for ml=find(marked_DistrStats&~isnan(distrStatsH'))
                set(distrStatsH(ml), 'LineWidth', linewidth_DistrStats(ml)+2);
            end
        end
    end



    function PointSelected(ObjectH, EventData,act_plot)
        
        C = get (get(ObjectH,'Parent'), 'CurrentPoint');
        xlim = get(get(ObjectH,'Parent'),'XLim');
        ylim = get(get(ObjectH,'Parent'),'YLim');
        Cnorm = [(C(1,1) - xlim(1))/(xlim(2)-xlim(1)) (C(1,2) - ylim(1))/(ylim(2)-ylim(1))];
        xnorm = (xdataCell{act_plot} - xlim(1))/(xlim(2)-xlim(1));
        ynorm = (ydataCell{act_plot} - ylim(1))/(ylim(2)-ylim(1));
        dst = sqrt((Cnorm(1)-xnorm).^2+(Cnorm(2)-ynorm).^2);
        
        %         dst = sqrt((C(1,1)-xdataCell{act_plot}).^2+(C(1,2)-ydataCell{act_plot}).^2);
        [~, idx] = min(dst);
        
        if isempty(marked_points)
            marked_points = zeros(max(cellfun(@(x) numel(x),xdataCell)),no_plots+1);
            markh = NaN(max(cellfun(@(x) numel(x),xdataCell)),no_plots+1);

        end
        
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
            
            markh(idx,act_plot) = plot(xdataCell{act_plot}(idx),ydataCell{act_plot}(idx),'s','Color',linecolor(act_plot,:),'MarkerSize',8); %get(ObjectH,'MarkerSize')+4
            uistack(markh(idx,act_plot),'bottom')
        end
    end

    function LineDelete(ObjectH, EventData,act_plot)
        lineH = update_LineHandle(figure_windows);
        xd = get(ObjectH,'XData');
        yd = get(ObjectH,'YData');
        
        xd2 = NaN; yd2 = NaN;
        if ~isempty(gco) && isprop(gco,'XData') && isprop(gco,'YData')
            xd2 = get(gco,'XData');
            yd2 = get(gco,'YData');
        end
        line_ud = get(ObjectH,'UserData');
        
        if isequal(xd,xd2) && isequal(yd,yd2)
            set(gco,'UserData',get(ObjectH,'UserData'));
            set(gco,'DeleteFcn',{@LineDelete,act_plot});
            lineH(act_plot) = gco;
            %             line_ud.listenerH1 =  addlistener(lineH(act_plot), 'LineWidth', 'PostSet', @(hAxes, eventData) setLineWidth(hAxes, eventData,act_plot,lineH(act_plot)));
            %             line_ud.listenerH2 =  addlistener(lineH(act_plot), 'Color', 'PostSet', @(hAxes, eventData) setColor(hAxes, eventData,act_plot,lineH(act_plot)));
            if isprop(gco,'Color') && isprop(ObjectH,'Color')
                set(gco,'Color',get(ObjectH,'Color'));
            elseif isprop(gco,'FaceColor') && isprop(ObjectH,'Color')
                set(gco,'FaceColor',get(ObjectH,'Color'));
            end
            set(ObjectH,'UserData',line_ud);
        else
            fig_ud = get(figure_windows,'UserData');
            fig_ud.act_set(fig_ud.act_set==line_ud.act_set) = [];
            set(figure_windows,'UserData',fig_ud);
            
            if ~isempty(line_ud) && isstruct(line_ud) && isfield(line_ud,'significanceMarkerH')
                if any(ishghandle(line_ud.significanceMarkerH))
                    delete(line_ud.significanceMarkerH(ishghandle(line_ud.significanceMarkerH)))
                end
            end
            if ~isempty(distrStatsH) && ishghandle(distrStatsH(act_plot))
                delete(distrStatsH(act_plot))
                distrStatsH(act_plot)=NaN;
            end
            if any(marked_points(:,act_plot)>0)
                delete(markh(ishandle(markh(:,act_plot)),act_plot));
                %             marked_points(:,act_plot) = [];
            end
        end
        
    end

    function  distrStatDelete(obj,event_obj)
        if ~any(ishandle(distrStatsH) & ishghandle(distrStatsH)) && ...
                numel(stat_menuH)>=6 && ~isempty(stat_menuH{6}) && ishghandle(stat_menuH{6})
            set(stat_menuH{6},'Enable','off');
        end
    end

    function remSignMarkers(obj,event_obj,sign_handle)
        if any(ishghandle(sign_handle))
            delete(sign_handle(ishghandle(sign_handle)));
        end
        if ~strcmpi(dm,'dynamicMapPolar') && ~strcmpi(dm,'MapPolar')
            ud = get(gca,'Userdata');
            if ~isempty(ud) && isstruct(ud)
                ud.ypos_signmarker = [];
                
                %                 ud = rmfield(ud,'ypos_signmarker');
                set(gca,'Userdata',ud);
            else
                set(gca,'Userdata',[]);
            end
            ylim1 = get(gca,'YLim');
            set(gca, 'ylimmode', 'auto')
            %         ylim2 = get(gca,'YLim');
            %         set(gca,'YLim',[ylim1(1) ylim2(2)]);
            %         ylim1 = get(gca,'YLim');
            %         set(gca,'YLim',[ylim1(1) ylim1(2)+(diff(ylim1)*0.05)]);
        end
    end

    function txt = datacursorupdatefcn(empt,event_obj)
        % Customizes text of data tips
        
        pos = get(event_obj,'Position');
        pval = get(get(event_obj,'Target'),'UserData');
        if ~isempty(pval) && isstruct(pval) && isfield(pval,'p')
            if isfield(pval,'test_type') && ~isempty(pval.test_type)
                if strcmpi(dm,'dynamicMapPolar') || strcmpi(dm,'MapPolar')
                    txt = {pval.test_type;['p = ',num2str(pval.p)]};
                else
                    txt = {pval.test_type;['p = ',num2str(pval.p(pos(1)==xticks))]};
                end
            else
                txt = ['p = ',num2str(pval.p(pos(1)==xticks))];
            end
        else
            txt = {['X: ',num2str(pos(1))],...
                ['Y: ',num2str(pos(2))]};
            
        end
    end

    function setLineWidth(hAxes, eventData,act_ln,lineHandle)
        if plotedit(figure_windows(end),'isactive')
            linewidth(act_ln) = get(lineHandle,'LineWidth');
        end
        
    end

    function setColor(hAxes, eventData,act_ln,lineHandle)
        ud = get(lineHandle,'UserData');
        if ~isempty(ud) && isfield(ud,'significanceMarkerH')
            c1 = get(ud.significanceMarkerH,'MarkerEdgeColor');
            c2 = get(ud.significanceMarkerH,'MarkerFaceColor');
            if isequal(linecolor(act_ln,:),c1)
                set(ud.significanceMarkerH,'MarkerEdgeColor',get(lineHandle,'Color'))
            elseif isequal(linecolor(act_ln,:),c2)
                set(ud.significanceMarkerH,'MarkerFaceColor',get(lineHandle,'Color'))
            end
        end
        if ~isempty(distrStatsH) && ishandle(distrStatsH(act_ln))
            set(distrStatsH(act_ln),'Color',get(lineHandle,'Color'))
        end
        linecolor(act_ln,:) = get(lineHandle,'Color');
    end

    function [src,stat_menuH] = addCallbacksToFigure(src)
        
        ud2 = get(src,'UserData');
        dm2 = ud2.display_mode;
        
        if isfield(ud2,'act_set') && ~isempty(ud2.act_set) && all(ud2.act_set(1)==ud2.act_set)
            act_set2 = floor(ud2.act_set(1));
        end
        
        stat_menuH = cell(1,16);
        

        if strcmpi(dm2,'MapPolar') || strcmpi(dm2,'dynamicMapPolar')
            stat_menuH{1} = uimenu(src,'Label','Statistics');
            stat_menuH{3} = uimenu(stat_menuH{1},'Label','Statistical Signifcance');
            
            stat_menuH{8} = uimenu(stat_menuH{3},'Label','Apply test','Enable','off');
            
            test_list = idSocial_auxiliaries_statisticalSignificance;
            for tl2 = 1:numel(test_list)
                spec_options = idSocial_auxiliaries_statisticalSignificance([],test_list{tl2});
                stat_menuH{4} = uimenu(stat_menuH{8},'Label',test_list{tl2},'Callback',{@setStatTest,spec_options,src});
            end
            
            stat_menuH{9} = uimenu(stat_menuH{1},'Label','Select data');
            stat_menuH{14} = uimenu(stat_menuH{1},'Label','Remove significance markers','Callback',{@remSignMarkers,sign_handle});
            
            ticks_actSet = input_data.(funcstring).results{8,act_set2+1};
            no_sets = size(input_data.(funcstring).results,2)-1;
            sign_sets = false(1,no_sets);
            for as = 1:no_sets
                if as~=act_set2
                    xt = input_data.(funcstring).results{8,as+1};
                    sign_sets(as) = isequal(ticks_actSet{1},xt{1}) & isequal(ticks_actSet{2},xt{2});
                end
            end
            
            stat_menuH{7} = NaN(1,no_sets);
            
            % If random controls exist:
            if size(input_data.(funcstring).results(:,act_set2+1),1)>9 && ~isempty(input_data.(funcstring).results{10,act_set2+1})
                stat_menuH{7}(no_sets+1) = uimenu(stat_menuH{9},'Label','Control','Callback',{@setSetSelect,no_sets+1});
            end
            
            for as = 1:no_sets
                if sign_sets(as)
                    stat_menuH{7}(as) = uimenu(stat_menuH{9},'Label',['Set ' num2str(as)],'Callback',{@setSetSelect,as});
                end
            end
            
            stat_menuH{13} = uimenu(stat_menuH{9},'Label','Clear selection','Callback',@clearSetSelection);
            
            stat_menuH{16} = datacursormode(src);
            set(stat_menuH{16},'UpdateFcn',@datacursorupdatefcn)
            
            
            
        else
            
            act_sets_ud2 = floor(ud2.act_set);
            distr_statisics_cell2 = cell(numel(act_sets_ud2),1);
            for k5 = 1:numel(act_sets_ud2)
                pm2 = input_data.(funcstring).(['Set' num2str(act_sets_ud2(k5))]).plot_mode;
                % For distributions:
                distr_statisics_cell2{k5} = pm2.statistics{1};
            end
            if numel(distr_statisics_cell2)>0 && all(cellfun(@(x) strcmp(distr_statisics_cell2{1},x),distr_statisics_cell2))
                distribution_flag = ~isempty(strfind(lower(distr_statisics_cell2{1}),'hist'));
                distr_statistics = strrep(lower(distr_statisics_cell2{1}), 'hist+','');
            end
            
            stat_menuH{1} = uimenu(src,'Label','Statistics');
            stat_menuH{2} = uimenu(stat_menuH{1},'Label','Show statistics');
            stat_menuH{3} = uimenu(stat_menuH{1},'Label','Statistical Signifcance');

            
            stat_menuH{4} = uimenu(stat_menuH{2},'Label','Standard deviation','Callback',@addSTD);
            stat_menuH{5} = uimenu(stat_menuH{2},'Label','SEM','Callback',@addSEM);
            if distribution_flag
                stat_menuH{6} = uimenu(stat_menuH{2},'Label',['Distribution ' distr_statistics],'Callback',@addDistrStats);
            end
            stat_menuH{7} = uimenu(stat_menuH{2},'Label','Remove statistics markers','Callback',@remError);
            
            
            stat_menuH{8} = uimenu(stat_menuH{3},'Label','Apply test');
            test_list = idSocial_auxiliaries_statisticalSignificance;
            stat_menuH{15} = cell(1,numel(test_list));
            for tl1 = 1:numel(test_list)
                spec_options = idSocial_auxiliaries_statisticalSignificance([],test_list{tl1});
                stat_menuH{15}{tl1} = uimenu(stat_menuH{8},'Label',test_list{tl1},'Callback',{@setStatTest,spec_options,src});
            end
            stat_menuH{9} = uimenu(stat_menuH{3},'Label','Select data');
            stat_menuH{10} = uimenu(stat_menuH{9},'Label','Select lines','Callback',@setLineSelect);
            stat_menuH{11} = uimenu(stat_menuH{9},'Label','Select data points','Callback',@setPointSelect);
            if distribution_flag
                stat_menuH{12} = uimenu(stat_menuH{9},'Label',['Select distribution ' distr_statistics],'Callback',@setDistrStatsSelect,'Enable','off');
            end
            stat_menuH{13} = uimenu(stat_menuH{9},'Label','Clear selection','Callback',@clearSelection);
            stat_menuH{14} = uimenu(stat_menuH{3},'Label','Remove significance markers','Callback',{@remSignMarkers,sign_handle});
            
                        
            
            stat_menuH{16} = datacursormode(src);
            set(stat_menuH{16},'UpdateFcn',@datacursorupdatefcn)
            
                    
            
            
        end
        ud2.stat_menuH = stat_menuH;
        set(src,'UserData',ud2);
        set(src,'CloseRequestFcn',{@closeFig})
    end

    function src=removeCallbacksToFigure(src)
%         mh = findobj(figure_windows,'Label','Statistics');
%         delete(mh);
        ud2 = get(src,'UserData');
        if ~isempty(ud2) && isfield(ud2,'stat_menuH')
            for k2 = numel(ud2.stat_menuH):-1:1
                if ~isempty(ud2.stat_menuH{k2})%% && ishghandle(ud2.stat_menuH{k2})
                    if iscell(ud2.stat_menuH{k2})
                        for sc = 1:numel(ud2.stat_menuH{k2})
                            if ~isempty(ud2.stat_menuH{k2}{sc})
                                delete(ud2.stat_menuH{k2}{sc});
                            end
                         end
                    else
                        if IsGraphicHandle(ud2.stat_menuH{k2})
                            delete(ud2.stat_menuH{k2})
                        end
                    end
                    ud2.stat_menuH{k2} = [];
                end
            end
        end
        
        % DeleteFcn of lines
        delH = findall(src,'-function',@(x) iscell(get(x,'DeleteFcn')));
        for dn = 1:numel(delH)
            df = get(delH(dn),'DeleteFcn');
            if iscell(df) && strcmp(func2str(df{1}),'idSocial_plotUI/LineDelete')
%                 delete(df)
                set(delH,'DeleteFcn','')
            end
        end
        
        menuH = findall(figure_windows,'Label','Statistics');
        delete(menuH);
        
        set(src,'CloseRequestFcn','')
        duh = datacursormode(src);
        set(duh,'UpdateFcn','')
        set(src,'UserData',ud2);
    end

    function [plotfh, input_data2, act_fig2, statsM] = newFigure(input_data2,funcstring,act_set)
        
        if nargin <3 || isempty(act_set)
            act_set = [];
        end
        
        
        plotfh = figure('Units','normalized','Position',[.45 .3 .5 .5],'NumberTitle','off','Color','w');
        
        fn2 = fieldnames(input_data2.(funcstring));
        maxFig2 = 0;
        for k2=1:size(fn2,1)
            if ~isempty(strfind(fn2{k2},'Fig'))
                maxFig2 = max(maxFig2,str2double(strrep(fn2{k2},'Fig','')));
            end
        end
        act_fig2 = maxFig2+1;
        
        [~, temp_savepath2]= idSocial_recursiveGetOptionsFromOptionsCell(input_data2.options,'temp_savepath');
        
        
        
        set(plotfh,'Name', ['Plot Results ' num2str(act_fig2)])
        ud.act_fig = act_fig2;
        ud.date = datestr(now, 30);
        ud.act_set = act_set;
        ud.save_path =  [temp_savepath2 funcstring '_figure' num2str(act_fig2) '.fig'];
        if ~isempty(act_set)
            pm2 = input_data2.(funcstring).(['Set' num2str(floor(act_set))]).plot_mode;
            ud.display_mode = pm2.display_mode;
        else
            ud.display_mode = [];
        end
        set(plotfh,'UserData',ud)
        
        input_data2.(funcstring).(['Fig' num2str(act_fig2)]).date = ud.date;
%         input_data2.(funcstring).(['Fig' num2str(act_fig2)]).handle = plotfh;
        input_data2.(funcstring).(['Fig' num2str(act_fig2)]).path = ud.save_path;
        
        if verLessThan('matlab','8.2')
            set(plotfh,'CloseRequestFcn','')
            plotfh = removeCallbacksToFigure(plotfh);
            saveas(plotfh, input_data2.(funcstring).(['Fig' num2str(act_fig2)]).path,'fig');
            set(plotfh,'CloseRequestFcn',{@closeFig})
%             plotfh = addCallbacksToFigure(plotfh);
            [plotfh,statsM] = addCallbacksToFigure(plotfh);
        else
            set(plotfh,'CloseRequestFcn','')
            plotfh = removeCallbacksToFigure(plotfh);
            savefig(plotfh, input_data2.(funcstring).(['Fig' num2str(act_fig2)]).path);
            set(plotfh,'CloseRequestFcn',{@closeFig})
%             plotfh = addCallbacksToFigure(plotfh);
            [plotfh,statsM] = addCallbacksToFigure(plotfh);
        end
        
        
    end
    function closeFig(src,evnt)
        
        uds2 = get(src,'UserData');
        save_path = uds2.save_path;
        if verLessThan('matlab','8.2')
            src = removeCallbacksToFigure(src);
            saveas(src, save_path,'fig');
        else
            src = removeCallbacksToFigure(src);
            pcell = findAndDeleteFunctionHandles(src);
            savefig(src, save_path);
%             dbstop if warning
        end
        delete(src)
        
    end
    function [pcell] = findAndDeleteFunctionHandles(fh)
        pnames = properties(fh);
        fcnidx  = cellfun(@(x) ~isempty(x),strfind(lower(pnames),'fcn')) | cellfun(@(x) ~isempty(x),strfind(lower(pnames),'callback'));
        
        pnamesFcn = pnames(fcnidx);
        pvals = get(fh,pnamesFcn)';
        % pvals = pvals(~isempty(pvals));
        
        handleidx  = cellfun(@(x) IsGraphicHandle(get(fh,x)),pnames,'UniformOutput',false);
        handleidx  = cellfun(@(x) any(x(:)),handleidx);
        
        hvals = get(fh,pnames(handleidx))';
        
        % p = ancestor(fh,'figure')
        pcell = [];
        try
        if ~all(cellfun(@(x) isempty(x),pvals));
            cellidx = cellfun(@(x) iscell(x),pvals);
            pvals(cellidx) = cellfun(@(x) x{1},pvals(cellidx),'UniformOutput',false);
            pnamesFcn = pnamesFcn(cellfun(@(x) ~isempty(x),pvals));

            pvals = pvals(cellfun(@(x) ~isempty(x),pvals));
            pstrings = cellfun(@(x) func2str(x),pvals,'UniformOutput',false);
            pcell = [pnamesFcn pstrings];
        end
        if ~isempty(pcell)
%             disp(pcell)
            set(fh,'WindowStyle','normal')
            % This should work in both HG1 and HG2:
            hManager = uigetmodemanager(fh);
            try
                set(hManager.WindowListenerHandles, 'Enable', 'off');  % HG1
            catch
                [hManager.WindowListenerHandles.Enabled] = deal(false);  % HG2
            end
            for k2 = 1:size(pcell,1)
                set(fh,pcell{k2,1},'')
            end
%             keyboard
            pvals2 = pvals(cellfun(@(x) ~isempty(x), pvals));
            %             delete(pvals(cellfun(@(x) ~isempty(x), pvals)))
        end
        
        catch
            keyboard
        end
        % delete(pvals(~isempty(pvals)))
    end

%     function Data = RemoveGraphicHandle(Data)
%         % Thanks to Bruno Luong, From https://www.mathworks.com/matlabcentral/newsreader/view_thread/342779
%         % Data = RemoveGraphicHandle(Data)
%         % Remove new graphics objects in a structure and/or cell
%         %
%         % Goal: Remove new graphic objects presented in Data, otherwise it crashs
%         % when loading the Data under earlier versions of MATLAB
%         % Note: Only use this function if NewGraphicsSystem() returns TRUE
%         
%         if ~NewGraphicsSystem()
%             return
%         end
%         
%         verbose = false;
%         
%         if isstruct(Data)
%             fn = fieldnames(Data);
%             for i = 1:numel(Data)
%                 for k = 1:length(fn)
%                     x = Data(i).(fn{k});
%                     if isstruct(x) || iscell(x) || IsNewGraphicHandle(x) || isa(x,'function_handle')
%                         Data(i).(fn{k}) = RemoveGraphicHandle(x);
%                     end
%                 end
%             end
%         elseif iscell(Data)
%             for k = 1:length(Data)
%                 x = Data{k};
%                 if isstruct(x) || iscell(x) || IsNewGraphicHandle(x) || isa(x,'function_handle')
%                     Data{k} = RemoveGraphicHandle(x);
%                 end
%             end
%         elseif isa(Data,'function_handle') && ...
%                 strcmp(getfield(functions(Data),'type'),'anonymous')
%             Data = func2str(Data);
%         elseif IsNewGraphicHandle(Data)
%             if verbose
%                 fprintf('Remove ');
%                 disp(Data)
%             end
%             % Replace with empty
%             Data = [];
%         end
%         
%     end % RemoveGraphicHandle
%
% %%
    

end

