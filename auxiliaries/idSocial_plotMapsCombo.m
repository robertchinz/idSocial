function subplot_handle=idSocial_plotMapsCombo(mapcell, optioncell,summary_titlestring,project_path)
% Plots the maps produces by idSocial_dynamics_makeMaps. 'mapno' determines
% the number of the map to be plotted, with mapno=1: heatmap, mapno=2...
% number of the map as defined by the index in 'optioncell{row_act,col_act}.mapstring'.
if nargin<3 || isempty(summary_titlestring)
    summary_titlestring='';
end

if iscell(summary_titlestring)
    summary_titlestring_complete=summary_titlestring;
    summary_titlestring=summary_titlestring{1};
else 
    summary_titlestring_complete=summary_titlestring;
end
    
%% Some design settings
columnlabelfontsize=14;
rowlabelfontsize=14;
titlefontsize=16;
axislabelsize=14;

%% Plots


uppercut=0; % This will cut off the upper and bottom x% of the z-values when calculating the color limits for the maps.


rows_subplot=size(mapcell,1);
cols_subplot=size(mapcell,2);
titleposition=[.045 .050 .85 .88];

% Plot heatmaps for all indices in mapcell
screen_size = get(0, 'ScreenSize');

% hf=figure('Position',[screen_size(3)-screen_size(4) 50  screen_size(4)-200 screen_size(4)-200]);

hf=figure;
set(gcf,'Tag',mfilename,'Name',summary_titlestring);

set(gcf,'renderer','OpenGL') % This is necessary to use colormaps>256, see http://www.mathworks.com/matlabcentral/answers/47458
set(gcf,'Units','Normal')
load('colormap4InteractionMaps.mat')
clims_all=NaN(rows_subplot,cols_subplot,2);
clims_accum=zeros(rows_subplot*cols_subplot+1,2);
pos_subplot=0;

cmap1 = jet(256);
cmap2 = colormap4InteractionMaps;
cmap_all=cat(1,cmap1,cmap2);

for row_act=1:size(mapcell,1)
    for col_act=1:size(mapcell,2)
        if size(mapcell{row_act,col_act},1)==1 && size(mapcell{row_act,col_act},2)==1 && ishandle(mapcell{row_act,col_act})
            hctr=findobj([mapcell{row_act,col_act}],'-regexp','Tag','meantimemap_contour_*');
            if ~isempty(hctr)
                parsetag=textscan(get(hctr,'Tag'),'%s','Delimiter','_');
                cmap3=jet(str2double(parsetag{1}{end})*2);
                cmap_all= cat(1,cmap_all,cmap3);
                break;
            end
        end
    end
end

% Define a colormap that consists of several separate colormaps.
% Idea from http://www.mathworks.es/support/solutions/en/data/1-GNRWEH/index.html
cmap_subplot=NaN(rows_subplot,cols_subplot);
for row_act=1:rows_subplot
    for col_act=1:cols_subplot
        pos_subplot=pos_subplot+1;
        map_act=mapcell{row_act,col_act};
        if min(map_act(:))<0
            cmap_subplot(row_act,col_act)=2;
            clim=max(abs(prctile(map_act(:),uppercut)),abs(prctile(map_act(:),100-uppercut)));
            
            clims_all(row_act,col_act,:)=[-clim clim];
            clims_accum(pos_subplot,:)=[clim 2*clim];
            
        else
            cmap_subplot(row_act,col_act)=1;
            clim=max(map_act(:));
            clims_all(row_act,col_act,:)=[0 clim];
            clims_accum(pos_subplot,:)=[0 clim];

        end
    end
end

pos_subplot=0;%cols_subplot;
subplot_handle=NaN(rows_subplot*cols_subplot,1);

for row_act=1:rows_subplot
    for col_act=1:cols_subplot
        
        
        
        pos_subplot=pos_subplot+1;
        if ~isempty(mapcell{row_act,col_act}) && ~all(ishandle(mapcell{row_act,col_act}(:)))
            map_act=mapcell{row_act,col_act};
            try
            if min(optioncell{row_act,col_act}.edges{1})<0
                noticsx=ceil(max(optioncell{row_act,col_act}.edges{1})*2);
            else
                noticsx=ceil(max(optioncell{row_act,col_act}.edges{1}));
            end
            catch
                keyboard
            end
            if min(optioncell{row_act,col_act}.edges{2})<0
                noticsy=ceil(max(optioncell{row_act,col_act}.edges{2})*2);
            else
                noticsy=ceil(max(optioncell{row_act,col_act}.edges{2}));
            end
            
            if  ~isempty(strfind(lower(optioncell{row_act,col_act}.heatmapstring{1}),'angle'))
                ticsx=round(optioncell{row_act,col_act}.edges{1}(1)):round(optioncell{row_act,col_act}.edges{1}(end));
                xtic=round((ticsx-optioncell{row_act,col_act}.edges{1}(1))*size(map_act,2)/(optioncell{row_act,col_act}.edges{1}(end)-optioncell{row_act,col_act}.edges{1}(1)));
                xtic(1)=1;
            else
                ticsx=optioncell{row_act,col_act}.edges{1}(round(optioncell{row_act,col_act}.edges{1})==optioncell{row_act,col_act}.edges{1});
                xtic=round(0:size(map_act,2)/noticsx:size(map_act,2));
                xtic(1)=1;
                
                %     xtic=xtic(2:end-1);
                %     ticsx=ticsx(2:end-1);
                
            end
            if  size(optioncell{row_act,col_act}.heatmapstring,1)==1 && ~isempty(strfind(lower(optioncell{row_act,col_act}.heatmapstring{1}),'angle')) || size(optioncell{row_act,col_act}.heatmapstring,1)==2 &&~isempty(strfind(lower(optioncell{row_act,col_act}.heatmapstring{2}),'angle'))
                ticsy=round(optioncell{row_act,col_act}.edges{2}(1)):round(optioncell{row_act,col_act}.edges{2}(end));
                ytic=round((ticsy-optioncell{row_act,col_act}.edges{2}(1))*size(map_act,1)/(optioncell{row_act,col_act}.edges{2}(end)-optioncell{row_act,col_act}.edges{2}(1)));
                ytic(1)=1;
            else
                ticsy=optioncell{row_act,col_act}.edges{2}(round(optioncell{row_act,col_act}.edges{2})==optioncell{row_act,col_act}.edges{2});
                ytic=round(0:size(map_act,1)/noticsy:size(map_act,1));
                ytic(1)=1;
                %     ytic=ytic(2:end-1);
                %     ticsy=ticsy(2:end-1);
            end
            
            
            
            
            subplot_handle(pos_subplot)=subplot(rows_subplot,cols_subplot,pos_subplot);
            p = get(subplot_handle(pos_subplot), 'pos');
            p(2)=p(2) - 0.05;
            p(3:4) = p(3:4) + 0.02;
            set(subplot_handle(pos_subplot), 'pos', p);
            
            
                imagesc(((map_act-clims_all(row_act,col_act,1))/(clims_all(row_act,col_act,2)-clims_all(row_act,col_act,1)))*255+1+(cmap_subplot(row_act,col_act)-1)*256+1);
                axis square tight
          
            
            colormap(cmap_all)
            hact=gca;
            set(gca,'YDir','normal','YTickLabelMode','Manual','YTickMode','Manual')
            %         if optioncell{row_act,col_act}.edges{2}(1)==optioncell{row_act,col_act}.edges{1}(1) && optioncell{row_act,col_act}.edges{2}(end)==optioncell{row_act,col_act}.edges{1}(end)
            
            %         end
            set(gca,'XTickLabel',[],'YTickLabel',[])
            ylabel({' '; optioncell{row_act,col_act}.plot_axislabelstring{2}},'fontsize',rowlabelfontsize);
            
            set(gca,'YTick',ytic,'YTickLabel',ticsy,'FontSize',axislabelsize)
            
            
            set(gca,'CLim', [0 max(cmap_subplot(:))*256])
            
            title(optioncell{row_act,col_act}.plot_titlestring{1},'fontsize',columnlabelfontsize)
            title(optioncell{row_act,col_act}.plot_titlestring{1},'fontsize',columnlabelfontsize) % I put this line here twice, because if not sometimes the title of the first subplot is missing. Sometimes it is there. Seems to work now, but I do not know what is going on.
            
            xlabel(optioncell{row_act,col_act}.plot_axislabelstring{1},'fontsize',rowlabelfontsize)
            set(gca,'XTick',xtic,'XTickLabel',ticsx,'FontSize',axislabelsize)
            
            %         if optioncell{row_act,col_act}.edges{2}(1)==optioncell{row_act,col_act}.edges{1}(1) && optioncell{row_act,col_act}.edges{2}(end)==optioncell{row_act,col_act}.edges{1}(end)
            axis square tight
            %         end
            
            ch=colorbar;
            labelstring=cellstr(num2str((clims_all(row_act,col_act,1):diff(clims_all(row_act,col_act,:))/4:clims_all(row_act,col_act,2))', '%.2f'));
            set(ch,'YLim', [(cmap_subplot(row_act,col_act)-1)*256+1 (cmap_subplot(row_act,col_act))*256+1],...
                'YTick',(cmap_subplot(row_act,col_act)-1)*256+1:(256-1)/4 :(cmap_subplot(row_act,col_act))*256+1,...
                'YTickLabel', labelstring)
            
        elseif ~isempty(mapcell{row_act,col_act}) && all(ishandle(mapcell{row_act,col_act}(:)))
            try
                
      
                
                subplot_handle(pos_subplot)=subplot(rows_subplot,cols_subplot,pos_subplot);
                if strcmpi(get(mapcell{row_act,col_act},'Type'),'figure')
                    allaxes = get(mapcell{row_act,col_act}, 'Child');           %# hAx = gca;
                    figpos=get(mapcell{row_act,col_act},'Position');
                    ax =findobj(allaxes,'Type','axes','-not','Tag','legend','-not','Tag','Colorbar','-not','Tag','scribeOverlay');
                elseif strcmpi(get(mapcell{row_act,col_act},'Type'),'axes') 
                    figpos=get(mapcell{row_act,col_act},'Position');
                    ax=mapcell{row_act,col_act}; 
                    allaxes=mapcell{row_act,col_act};
                end
                % Identify orientation of the axes
                [az,el] = view;
                
                axpos=get(ax,'Position');
%                 asp_ratio=axpos(3)/axpos(4);
                asp_ratio=figpos(3)/figpos(4);

                subpos=get(subplot_handle(pos_subplot),'position');
                %%%
                set(subplot_handle(pos_subplot),'PlotBoxAspectRatioMode','manual','PlotBoxAspectRatio',[asp_ratio 1 1]);
                view(subplot_handle(pos_subplot),[az,el])
                
                
                
                set(subplot_handle(pos_subplot),'Xlim',get(ax,'XLim'))
                set(subplot_handle(pos_subplot),'Ylim',get(ax,'YLim'))
                set(subplot_handle(pos_subplot),'Zlim',get(ax,'ZLim'))
                
                
                %%%
                delete(subplot_handle(pos_subplot))
                
                %             set(ax, 'Position',get(subplot_handle(pos_subplot)),'Position');
                
                %%% If handle is map
                hctr=findobj([mapcell{row_act,col_act}],'-regexp','Tag','meantimemap_contour_*');
%                 if ~isempty(hctr)
%                     keyboard
%                 end
                
                imh = findobj(ax,'Type','image');
                if ~isempty(imh)
                    map_act=get(imh,'CData');
                    if ~isempty(hctr)
                        cmap_subplot(row_act,col_act)=3;
                        parsetag=textscan(get(hctr,'Tag'),'%s','Delimiter','_');
                        cmap_depth=str2double(parsetag{1}{end});
                        clims=[0 cmap_depth];
                        map_act=((map_act-clims(1))/(clims(2)-clims(1)))*(cmap_depth*2)+2*256;
                    elseif min(map_act(:))<0
                        cmap_subplot(row_act,col_act)=2;
                        clim=max(abs(prctile(map_act(:),uppercut)),abs(prctile(map_act(:),100-uppercut)));
                        
                        clims=[-clim clim];
                        clims_accum(pos_subplot,:)=[clim 2*clim];
                        map_act=((map_act-clims(1))/(clims(2)-clims(1)))*255+1+(cmap_subplot(row_act,col_act)-1)*256+1;
                    elseif min(map_act(:))>0
                        cmap_subplot(row_act,col_act)=1;
                        clim=max(map_act(:));
                        clims=[0 clim];
                        clims_accum(pos_subplot,:)=[0 clim];
                        map_act=((map_act-clims(1))/(clims(2)-clims(1)))*255+1+(cmap_subplot(row_act,col_act)-1)*256+1;
                    end
%                     map_act=((map_act-clims(1))/(clims(2)-clims(1)))*255+1+(cmap_subplot(row_act,col_act)-1)*256+1;
                    
                    
                    colormap(cmap_all)
                     
%                     cb=colorbar;%colorbar('peer',ax);
%                     set(cb,'YLim', [(cmap_subplot(row_act,col_act)-1)*256+1 (cmap_subplot(row_act,col_act))*256+1]);
                end
                
                ch=copyobj(ax,hf);
                imh2 = findobj(ch,'Type','image');
                set(imh2,'CData',map_act);
                set(ch,'PlotBoxAspectRatioMode','manual','PlotBoxAspectRatio',[asp_ratio 1 1]);
                %                 set(ch,'CLim', [0 max(cmap_subplot(:))*256])
                set(ch,'CLim', [0 length(cmap_all)])
                
                subpos(2)=subpos(2)-(1-titleposition(4)-titleposition(2));
                
                set(ch,'position',subpos);
                set(ch,'Visible','On')
                if ~isempty(imh)
                    
                    %Find correct colormap
                    

                    cb=colorbar;%colorbar('peer',ax);
                    clims_cmap=[(cmap_subplot(row_act,col_act)-1)*256 min((cmap_subplot(row_act,col_act))*256+1,length(cmap_all))];
                    set(cb,'YLim', clims_cmap);

                    no_tics=5;
                    ytic=clims_cmap(1):diff(clims_cmap)/(no_tics):clims_cmap(2);
%                     yticlabel=clims(1):diff(clims)/(no_tics-1):clims(2);
                    yticlabel=cellstr(num2str((clims(1):diff(clims)/(no_tics):clims(2))', '%.1f'));
                    set(cb,'YTickMode','manual','YTickLabelMode','manual','YTick',ytic,'YTickLabel',yticlabel);
                    
                end
                

                
                %% Find and put legend
                legax =findobj(allaxes,'Type','axes','Tag','legend');
                if ~isempty(legax)
                    
                    loc=get(legax,'Location');
                    pos=get(legax,'Position');
                    newleg=legend(ch,'Location',loc);

                    copy_to_ax=newleg;
                    copy_from_ax=legax;
                    fnames=fieldnames(get(copy_from_ax));
                    
                    for fn=1:length(fnames)
                        if isprop(copy_to_ax,fnames{fn}) ...
                            && ~strcmp(fnames{fn},'Location') ...
                            && ~strcmp(fnames{fn},'Position') ...
                            && ~strcmp(fnames{fn},'OuterPosition') ...
                            && ~strcmp(fnames{fn},'ActivePositionProperty') ...
                            && ~strcmp(fnames{fn},'Parent')
                            try % try-catch to avoid 'read-only' properties'
                             set(copy_to_ax,fnames{fn},get(copy_from_ax,fnames{fn}))
                            catch
                            end
                        end
                    end
         
                 
                end
                
%                 keyboard

            catch err
                keyboard
                close(hf);
                msg = sprintf('%s', ...
                    'At least one figure has already been closed.'...
                    );
                error(['MATLAB:' mfilename ':InvalidFigure'], msg);
            end
        elseif isempty(mapcell{row_act,col_act})
            subplot_handle(pos_subplot)=subplot(rows_subplot,cols_subplot,pos_subplot);
            set(subplot_handle(pos_subplot),'Visible','off')
        else
                        
                
        end
    end
end

ax=axes('Units','Normal','Position',titleposition,'Visible','off');
set(get(ax,'Title'),'Visible','on')
title(summary_titlestring_complete,'fontsize',titlefontsize);

pause(.3)
jFrame = get(handle(hf),'JavaFrame');
jFrame.setMaximized(true)
pause(2)
%% Save a copy. As jpeg and fig.

idSocial_saveFigure(summary_titlestring,project_path);
% %% Save a copy. As jpeg and fig.
% thispath=mfilename('fullpath');
% breakidx=strfind(thispath,'\');
% savepath=thispath(1:breakidx(end-1));
% datestring=datestr(now, 30);
% saveas(hf, [savepath 'figures\' summary_titlestring '_' datestring],'fig');
% saveas(hf, [savepath 'figures\' summary_titlestring],'fig');
% set(gcf,'Color',[.95,.95,.95]) 
% set(hf, 'PaperPositionMode','auto')
% print(hf,[savepath 'figures\' summary_titlestring],'-djpeg','-r600');
% print(hf,[savepath 'figures\' summary_titlestring '_' datestring],'-djpeg','-r600');
% 
% disp(['Saved ' savepath 'figures\' summary_titlestring '_' datestring '.jpg and ~.fig' ]);