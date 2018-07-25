function input_data=idSocial_auxiliaries_signDiffMaps4groups(input_data,mapno,visible)


%% Some design settings
columnlabelfontsize=14;
rowlabelfontsize=14;
titlefontsize=16;
axislabelsize=14;
gr_tr_fontsize=12;

uppercut=0;%3; % This will cut off the upper and bottom x% of the z-values when calculating the color limits for the maps.
pos_uppercut=0;

if nargin<2 || isempty(visible)
    visible='Off';
end

signlevel=0.05;
statistics_function='ttest';
statistics_function='bootstrap';
statistics_function='ranksum';
field='meanmaps';
no_groups=size(input_data,1);
no_trials=size(input_data,2);

temp=arrayfun(@(x)x.data.(field){mapno},input_data,'uniformoutput',false);

meanmaps_alltrials_allfish=cell(no_groups,1);
C2=cell(size(temp{1,1},1),size(temp{1,1},2),no_groups);
for group=1:no_groups
    for trial=1:no_trials
        temp2=cat(3,temp{group,:});
        meanmaps_alltrials_allfish{group}=nanmean(temp2,3);
        temp2=permute(temp2,[3,1,2]);
        for row=1:size(temp{1,1},1)
            for col=1:size(temp{1,1},2)
                
                C2{row,col,group}=temp2(:,row,col);
            end
        end
    end
end
%%

mapcell=cell(no_groups+1,no_groups+1);
H=cell(no_groups+1,no_groups+1);
for group1=2:no_groups+1
    for group2=group1+1:no_groups+1
        mapcell{group1,group2}= meanmaps_alltrials_allfish{group1-1}-meanmaps_alltrials_allfish{group2-1};
        switch statistics_function
            case 'bootstrap'
            H{group1,group2}=arrayfun(@(x,y) strucdyn_bootstrap(x{1},y{1},100000),C2(:,:,group1-1),C2(:,:,group2-1));
            
            case 'ttest'
            H{group1,group2}=arrayfun(@(x,y)ttest2(x{1},y{1},signlevel,[],'unequal'),C2(:,:,group1-1),C2(:,:,group2-1));
            otherwise
                [p H{group1,group2}]=arrayfun(@(x,y)ranksum(x{1},y{1},'alpha',signlevel),C2(:,:,group1-1),C2(:,:,group2-1));
        end
        H{group2,group1}=H{group1,group2};
        mapcell{group2,group1}=-mapcell{group1,group2};
    end
    mapcell{1,group1}=meanmaps_alltrials_allfish{group1-1};
    mapcell{group1,1}=meanmaps_alltrials_allfish{group1-1};
end
%%
gl = opengl('data');
if strfind(gl.Vendor,'ATI')
    opengl('software')
end
rows_subplot=no_groups+1;
cols_subplot=no_groups+1;

edgesx=input_data(group,1).options.edges{1};
ticsx=floor(edgesx(1)):ceil(edgesx(end));
edgesy=input_data(group,1).options.edges{2};
ticsy=floor(edgesy(1)):ceil(edgesy(end));

% Plot heatmaps for all focals and neighbours
screen_size = get(0, 'ScreenSize');
fh=figure('Position',[screen_size(3)-screen_size(4) 50  screen_size(4)-200 screen_size(4)-200],'Visible',visible);
set(gcf,'Tag',mfilename,'Name',['Difference maps and significance for all groups: ' input_data(group,1).options.plot_titlestring{1}]);


load('colormap4InteractionMaps.mat') % Load colormap


clim=-inf;
temp=cat(3,meanmaps_alltrials_allfish{:});
clim_meanmaps= max(clim,max(abs(prctile(temp(:),uppercut)),abs(prctile(temp(:),100-uppercut))));

clim=-inf;
mapsdiff=mapcell(2:end,2:end);
temp=cat(3,mapsdiff{:});
clim= max(clim,max(abs(prctile(temp(:),uppercut)),abs(prctile(temp(:),100-uppercut))));
clim_diffmaps=clim;


clims_all=NaN(rows_subplot,cols_subplot,2);
clims_accum=zeros(rows_subplot*cols_subplot+1,2);
pos_subplot=0;

cmap1 = jet(256);
cmap2 = colormap4InteractionMaps;
cmap_all=cat(1,cmap1,cmap2);


% Define a colormap that consists of several separate colormaps.
% Idea from http://www.mathworks.es/support/solutions/en/data/1-GNRWEH/index.html
cmap_subplot=NaN(rows_subplot,cols_subplot);
for group1=1:rows_subplot
    for group2=1:cols_subplot
        pos_subplot=pos_subplot+1;
        map_act=mapcell{group1,group2};
        if min(map_act(:))<0
            cmap_subplot(group1,group2)=2;
            clim=clim_diffmaps;%max(abs(prctile(map_act(:),uppercut)),abs(prctile(map_act(:),100-uppercut)));
            
            clims_all(group1,group2,:)=[-clim clim];
            clims_accum(pos_subplot,:)=[clim 2*clim];
            
        else
            cmap_subplot(group1,group2)=1;
            clim=clim_meanmaps;%max(map_act(:));
            clims_all(group1,group2,:)=[0 clim];
            clims_accum(pos_subplot,:)=[0 clim];
            
        end
    end
end


pos_subplot=(no_groups+1)*(no_groups+1)+1;%cols_subplot;
subplot_handle=NaN((no_groups+1)*(no_groups+1),1);

for group1=no_groups+1:-1:1
    for group2=no_groups+1:-1:1
        pos_subplot=pos_subplot-1;
        
        
        %             group1=floor((pos_subplot-1)/cols_subplot)+1;
        %             group2=mod(pos_subplot-1,cols_subplot)+1;
        
        if ~isempty(mapcell{group1,group2})
            try
                subplot_handle(pos_subplot)=subplot(rows_subplot,cols_subplot,pos_subplot);
            catch
                keyboard
            end
            
            p = get(subplot_handle(pos_subplot), 'pos');
            p(2)=p(2) - 0.05;
            p(3:4) = p(3:4) + 0.02;
            set(subplot_handle(pos_subplot), 'pos', p);
            map_act=mapcell{group1,group2};
            if group1~=group2 && group2>1 && group1>1
                
                
                %%%%%%%
                
                %                             him= imagesc(edgesx,edgesy,mapcell{group1,group2},color_limits );
                him= imagesc(edgesx,edgesy,((map_act-clims_all(group1,group2,1))/(clims_all(group1,group2,2)-clims_all(group1,group2,1)))*255+1+(cmap_subplot(group1,group2)-1)*256+1);
                alphaval=.4;
                Halpha=double(H{group1,group2}); Halpha(Halpha==0)=alphaval;
                set(him,'AlphaDataMapping','none','AlphaData',Halpha)
                colormap(cmap_all)
                hact=gca;
                set(gca,'YDir','normal','YtickLabelMode','Manual','YtickMode','Manual')
                set(gca,'CLim', [0 max(cmap_subplot(:))*256])
                axis square tight
                set(gca,'XtickLabel',[],'YtickLabel',[])
                hold on
                
                xc=edgesx(1):(edgesx(end)-edgesx(1))/(size(H{group1,group2},2)-1):edgesx(end);
                yc=edgesy(1):(edgesy(end)-edgesy(1))/(size(H{group1,group2},1)-1):edgesy(end);
                contr=contourc(xc,yc,double(H{group1,group2}),[.5 .5]);
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
                
            elseif group1==1 && group2==1
                % The following plot will be invsible, left the whole stuff in case I want to put something else
                
                imagesc(edgesx,edgesy,((map_act-clims_all(group1,group2,1))/(clims_all(group1,group2,2)-clims_all(group1,group2,1)))*255+1+(cmap_subplot(group1,group2)-1)*256+1);
                colormap(cmap_all)
                set(h,'Visible','off')
                set(gca,'Visible','off')
                %             colormap(cmap_all)
                %             hact=gca;
                %             set(gca,'YDir','normal','YtickLabelMode','Manual','YtickMode','Manual')
                %             axis square tight
                %             set(gca,'XtickLabel',[],'YtickLabel',[])
                %
                %             ylabel({' '; input_data(group,1).options.plot_axislabelstring{2}},'fontsize',rowlabelfontsize);
                %             set(gca,'Xtick',ticsx,'Ytick',ticsy,'YtickLabel',ticsy,'FontSize',axislabelsize)
                %
                %             % The following lines put the title 'Mean' to the first figure.
                %             % Normally, the first of the following lines should do the
                %             % trick, but for some reason, the title appears only in about
                %             % 30% of the cases. Creating a second axes object with the same
                %             % title helps, even though I do not know why.
                %             ht=title(hact,'Mean','fontsize',columnlabelfontsize);
                %             B = axes('Position',get(hact,'Position'),'YLim',get(hact,'YLim'),'XLim',get(hact,'XLim'));
                %             set(B,'Xtick',ticsx,'Ytick',ticsy,'YtickLabel',[],'XtickLabel',[],'FontSize',axislabelsize)
                %             axis tight square
                %             title(B,'Mean','fontsize',columnlabelfontsize);
                %             uistack(gca,'bottom')
                set(gca,'Visible','off')
                
            elseif group2==1 && group1>1
                
                %                 imagesc(edgesx,edgesy,meanmaps_alltrials_allfish{group1-1},color_limits);
                imagesc(edgesx,edgesy,((map_act-clims_all(group1,group2,1))/(clims_all(group1,group2,2)-clims_all(group1,group2,1)))*255+1+(cmap_subplot(group1,group2)-1)*256+1);
                colormap(cmap_all)
                set(gca,'CLim', [0 max(cmap_subplot(:))*256])
                
                hact=gca;
                set(gca,'YDir','normal','YtickLabelMode','Manual','YtickMode','Manual')
                
                axis square tight
                set(gca,'XtickLabel',[],'YtickLabel',[])
                
                hl=ylabel({' '; input_data(group,1).options.plot_axislabelstring{2}},'fontsize',rowlabelfontsize);
                set(gca,'Xtick',ticsx,'Ytick',ticsy,'YtickLabel',ticsy,'FontSize',axislabelsize)
                
                
                B = axes('Position',get(hact,'Position'),'YLim',get(hact,'YLim'),'XLim',get(hact,'XLim'));
                set(B,'Xtick',ticsx,'Ytick',ticsy,'YtickLabel',[],'XtickLabel',[],'FontSize',axislabelsize)
                axis tight square
                %             title(B,'Mean','fontsize',columnlabelfontsize);
                uistack(gca,'bottom')
            elseif group2>1 && group1==1
                %             if group1==1 && group2==2
                %             keyboard
                %         end
                %                 imagesc(edgesx,edgesy,meanmaps_alltrials_allfish{group2-1},color_limits);
                imagesc(edgesx,edgesy,((map_act-clims_all(group1,group2,1))/(clims_all(group1,group2,2)-clims_all(group1,group2,1)))*255+1+(cmap_subplot(group1,group2)-1)*256+1);
                set(gca,'CLim', [0 max(cmap_subplot(:))*256])
                colormap(cmap_all)
                hact=gca;
                set(gca,'YDir','normal','YtickLabelMode','Manual','YtickMode','Manual')
                axis square tight
                set(gca,'XtickLabel',[],'YtickLabel',[])
                
                %             hl=ylabel({' '; input_data(group,1).options.plot_axislabelstring{2}},'fontsize',rowlabelfontsize);
                set(gca,'Xtick',ticsx,'Ytick',ticsy,'YtickLabel',[],'FontSize',axislabelsize)
                
                
                B = axes('Position',get(hact,'Position'),'YLim',get(hact,'YLim'),'XLim',get(hact,'XLim'));
                set(B,'Xtick',ticsx,'Ytick',ticsy,'YtickLabel',[],'XtickLabel',[],'FontSize',axislabelsize)
                
                axis tight square
                %             title(B,'Mean','fontsize',columnlabelfontsize);
                uistack(gca,'bottom')
            else % Placeholder
                h=imagesc(edgesx,edgesy,((map_act-clims_all(1,2,1))/(clims_all(1,2,2)-clims_all(1,2,1)))*255+1+(cmap_subplot(1,2)-1)*256+1);
                set(gca,'CLim', [0 max(cmap_subplot(:))*256])
                colormap(cmap_all)
                set(h,'Visible','off')
                set(gca,'Visible','off')
                set(get(gca,'Title'),'Visible','on')
                set(get(gca,'YLabel'),'Visible','On')
                
                
            end
            
            
            %%% Put title, labels etc. depending on rows and columns
            
            %%%% Titles 1st row and 1st column, 2nd row
            if group1==1 && group2>1 %|| group1==2 && group2==2
                title(['Group ' num2str(group2-1)],'fontsize',columnlabelfontsize)
            elseif group1==2 && group2==1
                ht=title('Maps','fontsize',columnlabelfontsize);
                
            end
            
            %%%% Y-ticks for first column
            if group2==1 % && group1>1
                hylabel=ylabel( input_data(group,1).options.plot_axislabelstring{2},'fontsize',rowlabelfontsize);
                set(gca,'Xtick',ticsx,'Ytick',ticsy,'YtickLabel',ticsy,'FontSize',axislabelsize)
                
            end
            
            %%%% Y-axis label  for last column, first row
            
            if group1==1 && group2==cols_subplot && group1<rows_subplot
                
                ylabel(gca,'Maps','fontsize',rowlabelfontsize)
                set(gca,'yaxislocation','right','Ytick',[]);
                
                if group2==cols_subplot-1 && cols_subplot==2 % If N=2 put y-label on both sides
                    B = axes('Position',get(hact,'Position'),'YLim',get(hact,'YLim'),'XLim',get(hact,'XLim'));
                    set(B,'yaxislocation','left','Xtick',ticsx,'Ytick',ticsy,'YtickLabel',ticsy,'FontSize',axislabelsize)
                    axis tight square
                    ylabel(B, input_data(group,1).options.plot_axislabelstring{2},'fontsize',rowlabelfontsize)
                    uistack(gca,'bottom')
                end
                
            end
            
            %%%% Y-axis label  for last column
            
            if group1>1 && group2==cols_subplot && group1<rows_subplot || group1>1 && group2==cols_subplot-1 && group1==rows_subplot
                
                ylabel(gca,['Group ' num2str(group1-1)],'fontsize',rowlabelfontsize)
                set(gca,'yaxislocation','right','Ytick',[]);
                
                if group2==cols_subplot-1 && cols_subplot==2 % If N=2 put y-label on both sides
                    B = axes('Position',get(hact,'Position'),'YLim',get(hact,'YLim'),'XLim',get(hact,'XLim'));
                    set(B,'yaxislocation','left','Xtick',ticsx,'Ytick',ticsy,'YtickLabel',ticsy,'FontSize',axislabelsize)
                    axis tight square
                    ylabel(B, input_data(group,1).options.plot_axislabelstring{2},'fontsize',rowlabelfontsize)
                    uistack(gca,'bottom')
                end
                
            end
            
            %%%% X-axis label last row and last subplot in second last row
            if group1==rows_subplot || group1==rows_subplot-1 && group2==cols_subplot
                xlabel(input_data(group,1).options.plot_axislabelstring{1},'fontsize',rowlabelfontsize)
                set(gca,'Xtick',ticsx,'Ytick',ticsy,'XtickLabel',ticsx,'FontSize',axislabelsize)
                
            end
            
            if group1==1 && rows_subplot==2 && group2==2 % If N=2 put y-label on both sides
                B = axes('Position',get(hact,'Position'),'YLim',get(hact,'YLim'),'XLim',get(hact,'XLim'));
                set(B,'Xtick',ticsx,'Ytick',ticsy,'XtickLabel',ticsx,'FontSize',axislabelsize)
                axis tight square
                xlabel(B, input_data(group,1).options.plot_axislabelstring{1},'fontsize',rowlabelfontsize)
                uistack(gca,'bottom')
            end
            
            
            
            
            
        end
        %%%% Colorbars
        
        if group1==1 && group2==1
            subplot_handle(pos_subplot)=subplot(rows_subplot,cols_subplot,pos_subplot);
            map_act=mapcell{1,2};
            h=imagesc(edgesx,edgesy,((map_act-clims_all(1,2,1))/(clims_all(1,2,2)-clims_all(1,2,1)))*255+1+(cmap_subplot(1,2)-1)*256+1);
            set(gca,'CLim', [0 max(cmap_subplot(:))*256])
            
            set(h,'Visible','off')
            set(gca,'Visible','off')
            
            colormap(cmap_all)
            
            hc=colorbar;
            set(hc,'location','southoutside');
            labelstring=cellstr(num2str((clims_all(1,2,1):diff(clims_all(1,2,:))/4:clims_all(1,2,2))', '%.2f'));
            set(hc,'XLim', [(cmap_subplot(1,2)-1)*256+1 (cmap_subplot(1,2))*256+1],...
                'Xtick',(cmap_subplot(1,2)-1)*256+1:(256-1)/4 :(cmap_subplot(1,2))*256+1,...
                'XtickLabel', labelstring);
            title(hc,{'Means'},'fontsize',rowlabelfontsize)
            hc_pos=get(hc,'Position');
            hc_pos(1)=hc_pos(1)+.01;
            hc_pos(3:4)=hc_pos(3:4);%*1.1;
            set(hc,'Position',hc_pos)
        end
        
        if group1==rows_subplot && group2==cols_subplot
            subplot_handle(pos_subplot)=subplot(rows_subplot,cols_subplot,pos_subplot);
            map_act=mapcell{no_groups-1,no_groups};
            h=imagesc(edgesx,edgesy,((map_act-clims_all(no_groups-1,no_groups,1))/(clims_all(no_groups-1,no_groups,2)-clims_all(no_groups-1,no_groups,1)))*255+1+(cmap_subplot(no_groups-1,no_groups)-1)*256+1);
            set(gca,'CLim', [0 max(cmap_subplot(:))*256])
            
            set(h,'Visible','off')
            set(gca,'Visible','off')
            
            colormap(cmap_all)
            
            hc=colorbar;
            set(hc,'location','southoutside');
            labelstring=cellstr(num2str((clims_all(group1-1,group2,1):diff(clims_all(group1-1,group2,:))/4:clims_all(group1-1,group2,2))', '%.2f'));
            set(hc,'XLim', [(cmap_subplot(group1-1,group2)-1)*256+1 (cmap_subplot(group1-1,group2))*256+1],...
                'Xtick',(cmap_subplot(group1-1,group2)-1)*256+1:(256-1)/4 :(cmap_subplot(group1-1,group2))*256+1,...
                'XtickLabel', labelstring);
            title(hc,{'Differences'},'fontsize',rowlabelfontsize)
            hc_pos=get(hc,'Position');
            hc_pos(1)=hc_pos(1)+.01;
            hc_pos(3:4)=hc_pos(3:4);%*1.1;
            set(hc,'Position',hc_pos)
        end
    end
end
%         xpostxt=get(hylabel,'Position'); xpostxt=xpostxt(1);
%         ypostxt=get(htitle,'Position'); ypostxt=ypostxt(2);
%         text(xpostxt,ypostxt,['(' num2str(group) ',:)'],'fontsize',gr_tr_fontsize);
ax=axes('Units','Normal','Position',[-.33 .081 .85 .85],'Visible','off');
set(get(ax,'Title'),'Visible','on','horizontalAlignment', 'center')
title(['(:,:)'],'fontsize',gr_tr_fontsize);
ax=axes('Units','Normal','Position',[.095 .075 .85 .85],'Visible','off');
set(get(ax,'Title'),'Visible','on','horizontalAlignment', 'center')
htitle=title(input_data(group,1).options.plot_titlestring{1},'fontsize',titlefontsize);
annotation('textbox', [.7, .05, 1, 0], 'string', ['Statistical test: ' statistics_function ' at significance level ' num2str(signlevel)]);

%% Put color background

%% Put everything in structure
% Generate string for fieldname
filename=mfilename;
regstr=regexprep(input_data(group,1).options.plot_titlestring{1},'[^a-zA-Z0-9]','_');
scnstr=textscan(regstr,'%s','Delimiter','_');
titstr=[];

for p=1:size(scnstr{1},1)
    if ~isempty(scnstr{1}{p})
        scnstr{1}{p}(1)=upper(scnstr{1}{p}(1));
        titstr=[titstr scnstr{1}{p}(1:min(4,length(scnstr{1}{p})))];
    end
end


input_data(group,1).figures.([filename(10:end) '_' titstr '_spl']).handle=subplot_handle;
input_data(group,1).figures.([filename(10:end) '_' titstr]).handle=fh;

%% Save a copy. As jpeg and fig.

idSocial_saveFigure([filename(10:end) '_' titstr '_' num2str(group)],input_data(group,1).options.project_path)

input_data(group,1).figures.([filename(10:end) '_' titstr '_spl']).path=[filename(10:end) '_' titstr '_spl_' num2str(group)];
input_data(group,1).figures.([filename(10:end) '_' titstr]).path=[filename(10:end) '_' titstr '_' num2str(group) ];
%     end


