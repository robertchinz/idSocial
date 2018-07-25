function idSocial_plotGraph(input_data,act_method,plot_mode,plot_options)
% input_data.(act_method).output_plot.data is a N-D-array
% containing a M-by-M matrix (besides of singleton
% dimensions) with elements (i,j).
% Direction of arrows: (i,j)<0 (and (j,i)>0): Arrow from i
% to j.
% Vertical position of node i: Determined by
% Mean([-(i,1),-(i,2),...,-(i,M)])

%% Check input
if nargin<3 || ~isfield(plot_mode,'visible')|| isempty( plot_mode.visible)
    visible='On';
else
    visible=plot_mode.visible;
end
if nargin<3 || ~isfield(plot_mode,'display_mode') || isempty( plot_mode.display_mode)
    display_mode='plot2d';
else
    display_mode=input_data(1,1).(act_method).output_plot.display_mode;
end

if nargin<4 || isempty(plot_options)
    plot_options = [];
end

set(0,'DefaultPatchLineSmoothing','On')
set(0,'DefaultLineLineSmoothing','On')

%%
% Set options
options=input_data(1,1).(act_method).options;
% act_method=         options.act_method;
edges=              options.edges;
titlefontsize=      options.plot_titlefontsize;
tickfontsize=       options.plot_tickfontsize;
axislabelsize=      options.plot_axislabelsize;
labelfontsize=      options.plot_labelfontsize;
gr_tr_fontsize=     options.plot_gr_tr_fontsize;
axis_scaling=       options.plot_axis_scaling;
colororder=               options.plot_colororder;
linewidth=          options.plot_linewidth;
timeintervals_in_min=options.timeintervals_in_min;
axislabelstring=    options.plot_axislabelstring;
titlestring=        options.plot_titlestring;
info=               input_data(1,1).info;
framerate=          info.framerate;
plot_xticklabel=        options.plot_xticklabel;
plot_bootstrap_repetitions= options.plot_bootstrap_repetitions;
% plot_random=          options.plot_random;
spacing=    options.spacing;
significance_between_groups= options.significance_between_groups;
statistical_test_type=      options.statistical_test_type;
statistical_test_significance=                  options.statistical_test_significance;
linespec=       {'--'; '-.'; '-'};
group_names=    info.group_name;

%% Some standard parameters
output_plot=input_data(1,1).(act_method).output_plot;

no_groups=size(input_data,1);
no_trials=size(input_data,2);
max_no_fish=max(max(info.no_focals));
no_frames=  info.no_frames;
duration=           info.duration;
if isempty(timeintervals_in_min) || timeintervals_in_min> min(floor(min(info.no_frames./info.framerate)/60));
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


def_height_txt=.2;

%% Plot the stuff



timeidx=find(strcmpi(output_plot.dim_names_original,'time'));
groupidx=find(strcmpi(output_plot.dim_names_original,'group'));
edgesXidx=find(strcmpi(output_plot.dim_names_original,'edgeX'));
edgesYidx=find(strcmpi(output_plot.dim_names_original,'edgeY'));

%% Calculate y- (or color-) limits
lim=[0 0];


no_subplots=size(output_plot.data,1);
no_legends=size(output_plot.data,2);

v=[];
v_dev_up=[];
v_dev_down=[];
conf_intervall=cell(no_subplots,no_legends,2);
for sp=1:no_subplots
    for lg=1:no_legends
        v=[v output_plot.data{sp,lg}(:)'];
        if ~isempty(output_plot.data_sign) && ~isempty(output_plot.data_sign{sp,lg})
            switch output_plot.statistics_type{sp,lg}
                case {'Mean','mean','MEAN'}
                    v_dev_up=[v_dev_up output_plot.data_dev{sp,lg}(:)'./sqrt(output_plot.data_no_datapoints{sp,lg}(:)')];
                    v_dev_down=v_dev_up;
                case {'Median','median','MEDIAN'}
                    S.type='()';
                    S.subs=repmat({':'},[size(output_plot.cell_size,2),1]);
                    S.subs{output_plot.statistics_on_idx(sp,lg)}=1;
                    conf_intervall{sp,lg,1}=squeeze(subsref(output_plot.data_dev{sp,lg},S));
                    v_dev_up=[v_dev_up conf_intervall{sp,lg,1}(:)'];
                    
                    S.subs=repmat({':'},[size(output_plot.cell_size,2),1]);
                    S.subs{output_plot.statistics_on_idx(sp,lg)}=2;
                    conf_intervall{sp,lg,2}=subsref(output_plot.data_dev{sp,lg},S);
                    v_dev_down=[v_dev_down conf_intervall{sp,lg,2}(:)'];
            end
        else
            v_dev_up=0;
            v_dev_down=0;
        end
        
    end
    
end
maxv=max(v(:)+v_dev_up(:));
minv=min(v(:)-v_dev_down(:));

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
        
        lim=[0 maxv*1.1];
end
if lim(1)==lim(2); lim(1)=lim(1)-.1;  lim(2)=lim(2)+.1; end
%% Set xaxis

if output_plot.xaxis_idx==timeidx
    axislabelstring{1}='Time [min]';
end


%%


sprows=round(sqrt(no_legends));
spcols=ceil(sqrt(no_legends));

subplotstring={''};
for sp=1:max(1,no_subplots)
    
    if ~isempty(output_plot.subplot_idx)
        sps_act=output_plot.subplot_idx(sp,:);
        
        
        for ss=1:length(output_plot.subplot_idx{sp,1})
            if output_plot.subplot_idx{sp,1}(ss)==timeidx
                output_plot.subplotstring{sp}=['t= ' num2str(timevec(output_plot.subplot_idx{sp,2}(ss),1),'%.1f') '-' num2str(timevec(min(output_plot.subplot_idx{sp,2}(ss),length(timevec)),2),'%.1f') 'min'];
            elseif output_plot.subplot_idx{sp,1}(ss)==groupidx
                output_plot.subplotstring{sp}=group_names{sps_act{ss,2}};
            end
        end
        figname=[act_method ', ' output_plot.subplotstring{sp}];
    else
        figname=act_method;
    end
    
    
    titlestr_prefix=[];
    %     fh_raw=figure('Visible','off','units','normalized','outerposition',[.2 .1 .7 .8],'Color',[.95,.95,.95]);
    %     set(fh_raw,'Tag',[act_method '_raw'])
    %     set(fh_raw,'Name',figname);
    fh_raw = gcf;
    set(fh_raw,'Visible','on','units','normalized','outerposition',[.1 .1 .4 .8]);
    set(fh_raw,'Tag',[act_method '_raw'])
    set(fh_raw,'Name',figname);
    axis off
    no_subplots=size(output_plot.data,1);
    no_legends=size(output_plot.data,2);

    sprows=round(sqrt(no_legends));
    spcols=ceil(sqrt(no_legends));
    annotstr='Figure shows: ';
    sphandle=NaN(no_subplots,1);
    for lg=1:no_legends
        row_act=sprows-floor((lg-1)/spcols+1)+1;
        col_act=mod((lg-1),spcols)+1;
        
        sphandle(lg)=axes('Units','normalized','OuterPosition',[(col_act-1)*1/spcols def_height_txt+(row_act-1)*(1-def_height_txt)/sprows 1/spcols (1-def_height_txt)/sprows*.85]);
        
        %         subplot(sprows,spcols,lg)
        he=idSocial_delay2grafo(-squeeze(output_plot.data{sp,lg}),~isnan(squeeze(squeeze(output_plot.data{sp,lg}))),1,[],colororder);
%         axis off equal
        xlabel(axislabelstring{1},'fontsize',labelfontsize,'Color','k')
        ylabel(axislabelstring{2},'fontsize',labelfontsize,'Color','k');
        try
            
            if ~isempty(plot_options) && isfield(plot_options,'plot_titlestring')
                title(plot_options.plot_titlestring,'fontsize',titlefontsize,'Color','k');
            elseif ~isempty(output_plot.subplotstring)
                title({[titlestr_prefix titlestring{1} ' ' output_plot.subplotstring{sp}]; ...
                    output_plot.legendstring{sp,lg}},'fontsize',titlefontsize,'Color','k');
            elseif (isempty(plot_options) || isfield(plot_options,'plot_titlestring')) && isempty(output_plot.subplotstring)
                title({[titlestr_prefix titlestring{1} ];...
                    output_plot.legendstring{sp,lg}},'fontsize',titlefontsize,'Color','k');
            end
        catch
            keyboard
        end
        set(sphandle(lg),'UserData','top')
        
    end
    newfig=fh_raw;%  myaa;
    %       close(fh_raw)
    set(newfig,'Visible',visible, ...
        'Tag',[act_method '_raw'],'Name',figname);
    
    
    % Annotations and legend
    if max(1,no_subplots)>1
        annotstr=[annotstr 'Fig. ' num2str(sp) ') '];
    end
    for lg=1:no_legends
        annotstr=[annotstr output_plot.data_string{sp,lg} '; '];
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
    
    
end
% idSocial_saveFigure([act_method '_' fname],project_path)
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
        
        if new_height>0
            set(ax1(k),'OuterPosition',[posx new_posy width new_height])
        else
            new_height=act_pos(4);
            set(ax1(k),'OuterPosition',[posx new_posy width new_height])
        end
    end
end
end
