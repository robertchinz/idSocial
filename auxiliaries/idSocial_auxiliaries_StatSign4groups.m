function [H p CI stat sampstat]= ...
    idSocial_auxiliaries_StatSign4groups(val,stattest,statistical_test_significance,ax,colororder,bootstrap_repetitions,display_mode,bootstrap_statfunc,tail,compare_rand_only,mode_edges,mode_bandwidth,mode_method,sameSampleSize)
% Fomat val: (groups,trials,parts)
no_groups=size(val,1);

no_parts=size(val,3);

if nargin<7 || isempty(display_mode)
    display_mode='plot2d';
end
if nargin<6 || isempty(bootstrap_repetitions)
    bootstrap_repetitions=10000;
end

if nargin<8 || isempty(bootstrap_statfunc)
    bootstrap_statfunc='Mean';
end


if nargin<9 || isempty(tail)
    tail='both';
end

if nargin<10 || isempty(compare_rand_only)
    compare_rand_only=false;
end

if nargin<11 || isempty(mode_edges) 
    mode_edges = [];
end

if nargin<12 || isempty(mode_bandwidth)
    mode_bandwidth = []; 
end

if nargin<13 || isempty(mode_method)
    mode_method = []; 
end

if nargin<14 || isempty(sameSampleSize)
    sameSampleSize = true; 
end


if nargin<5 || isempty(colororder)
    colororder=[1 0 0; 0 1 0; 0 0 1; 1 1 0; 0 1 1; 1 0 1; .5 0 1; 1 0 0; 0 1 0; 0 0 1; 1 0 1; 1 .5 .5; 0 1 1; 1 1 0;  .25 .5 .5; 1 .5 0; 0 .6 0; .5 .2 1; .5 1 1; 1 1 0];
end


% if iscell(bootstrap_statfunc) && size(bootstrap_statfunc,2) == 2
%     bootstrap_statfunc_param = bootstrap_statfunc{1,2};
%     bootstrap_statfunc = bootstrap_statfunc{1,1};  
% end
%% Calculate statistical significance between groups

CI=[];
H=NaN(no_parts,no_groups,no_groups);
p=NaN(no_parts,no_groups,no_groups);
stat=NaN(no_parts,no_groups,no_groups,bootstrap_repetitions);
sampstat=NaN(no_parts,no_groups,no_groups);
if compare_rand_only
    group_index=[(1:no_groups/2)' (no_groups/2+1:no_groups)'];
else
    group_index=vectorsize2indexcombinations([no_groups no_groups]);
end 

for act_idx=1:size(group_index,1)
    %     for group2=group_index
    group=group_index(act_idx,1);
    group2=group_index(act_idx,2);
    if group~=group2
        switch stattest
            case 'ttest'
                [H(:,group,group2) p(:,group,group2)]=ttest2(squeeze(val(group,:,:)),squeeze(val(group2,:,:)),statistical_test_significance,tail,'unequal');
            case 'bootstrap'
                for step=1:no_parts
                    [H(step,group,group2) p(step,group,group2) CI stat(step,group,group2,:) sampstat(step,group,group2)]=idSocial_bootstrap(squeeze(val(group,:,step))',squeeze(val(group2,:,step))',bootstrap_repetitions,statistical_test_significance,bootstrap_statfunc, tail);
                end
            case 'polaviejaStat'
               
                for step=1:no_parts
%                     if isequal(squeeze(val(group,:,step)),squeeze(val(group2,:,step)))
%                         keyboard
%                     end
%                     if group==4 && group2==9; keyboard; end
                    try
                        [H(step,group,group2) p(step,group,group2) CI stat(step,group,group2,:) sampstat(step,group,group2)]=idSocial_polaviejaStat(squeeze(val(group,:,step))',squeeze(val(group2,:,step))',bootstrap_repetitions,statistical_test_significance,...
                            bootstrap_statfunc, tail, mode_edges, mode_bandwidth, mode_method,sameSampleSize);
                        
                    catch
                        keyboard
                    end
                end
            case 'polaviejaStatRandControls'
                
                for step=1:no_parts
                    %                     if isequal(squeeze(val(group,:,step)),squeeze(val(group2,:,step)))
                    %                         keyboard
                    %                     end
                    %                     if group==4 && group2==9; keyboard; end
                    try
                        [H(step,group,group2) p(step,group,group2) CI stat(step,group,group2,:) sampstat(step,group,group2)]=idSocial_polaviejaStat4RandomizedControls(squeeze(val(group,:,step))',squeeze(val(group2,:,step))',bootstrap_repetitions,statistical_test_significance,...
                            bootstrap_statfunc, tail, mode_edges, mode_bandwidth, mode_method);
                        
                    catch
                        keyboard
                    end
                end
            case 'kstest2'
                for step=1:no_parts
                    [p(step,group,group2) H(step,group,group2)]=kstest2(squeeze(val(group,~isnan(val(group,:,step)),step))',squeeze(val(group2,~isnan(val(group2,:,step)),step)),statistical_test_significance);
                end
            otherwise
                stattest='ranksum';
                vr = version('-release');
                vr_yy = str2double(vr(1:4));
                vr_aa = vr(end);
                if (vr_yy == 2012 && vr_aa > 'a') || vr_yy > 2012
                    for step=1:no_parts
                        if ~(all(isnan(val(group,:,step)))) && ~(all(isnan(val(group2,:,step))))
%                              keyboard
                            [p(step,group,group2) H(step,group,group2)]=ranksum(squeeze(val(group,~isnan(val(group,:,step)),step))',squeeze(val(group2,~isnan(val(group2,:,step)),step)),'alpha',statistical_test_significance,'tail',tail);
                        end
                    end
                else
                    warning([mfilename ': No tail options available for ranksum in Matlab versions < 2012b. Using both sided test!'])
                    disp([mfilename ': No tail options available for ranksum in Matlab versions < 2012b. Using both sided test!'])
                  
                    for step=1:no_parts
                        if ~(all(isnan(val(group,:,step)))) && ~(all(isnan(val(group2,:,step))))
                            [p(step,group,group2) H(step,group,group2)]=ranksum(squeeze(val(group,~isnan(val(group,:,step)),step))',squeeze(val(group2,~isnan(val(group2,:,step)),step)),'alpha',statistical_test_significance);
                        end
                    end
                end
                
        end
        
        %     end
    end
end

if nargin>3 && ~isempty(ax) && ishandle(ax)
    
    val_mean=squeeze(nanmean(val,2));
    val_std=squeeze(nanstd(val,[],2));
    
    yl=get(ax,'YLim');
    
    
    ax_yscale=diff(get(ax,'YLim'));
  
    
    if strcmp(display_mode,'hist')
        
        lineObj=findobj(get(ax,'Children'),'Type','line','-or','Type','stair');
        ydat=get(lineObj,'YData');
        if size(ydat,1)==1 && ~iscell(ydat)
            ydat={ydat};
        end
        ydat=ydat(cellfun(@(x) size(x,2)>1,ydat));
        pvaluefontsize=10;
        linewidth=2;
        for part=1:size(H,1)
            for group1=1:no_groups
                for group2=group1+1:no_groups
                    max_val_at_t=max(ydat{group1},[],2)';
                    yrange=yl(2)*.9-max_val_at_t;
%                     ypos_signmarker=max_val_at_t;
                    
                    
                    ypos_signmarker=max_val_at_t+nansum(nansum(nansum(squeeze(H(1:part,1:group1,1:group2)),1),2))/(nansum(H(:))/2)*yrange;
                    ypos_signmarker(H(part,group1,group2)==0 | isnan(H(part,group1,group2)))=NaN;
                    text((nanmean(val(group1,:,part))+nanmean(val(group2,:,part)))/2,ypos_signmarker,['p=' num2str(p(part,group1,group2),'%.2g')],...
                        'Color',colororder(group2,:),'HorizontalAlignment','center','VerticalAlignment','bottom','fontsize',pvaluefontsize,...
                        'BackgroundColor','w','Margin',1)
                    
                    plot([nanmean(val(group1,:,part)) nanmean(val(group2,:,part))],[ypos_signmarker ypos_signmarker],'-','Color',colororder(group1,:),'MarkerFaceColor',colororder(group2,:),'MarkerSize',8,'LineWidth',linewidth)
                    plot([nanmean(val(group1,:,part)) nanmean(val(group1,:,part))],yl,'-','Color',colororder(group1,:),'MarkerFaceColor',colororder(group1,:),'MarkerSize',8,'LineWidth',linewidth)
                    plot([nanmean(val(group2,:,part)) nanmean(val(group2,:,part))],yl,'-','Color',colororder(group2,:),'MarkerFaceColor',colororder(group2,:),'MarkerSize',8,'LineWidth',linewidth)
                   
                end
            end
        end
        
        
    else
        
        lineObj=findobj(get(ax,'Children'),'Type','line');
        if isempty(lineObj)
         lineObj=findobj(get(ax,'Children'),'Type','errorbar');
        end
        if isfield(get(lineObj(1)),'UData')
             ydat=cellfun(@(x,y) x+y,get(lineObj,'YData'),get(lineObj,'UData'),'UniformOutput',false);
        else
            ydat=get(lineObj,'YData');
        end
        if size(ydat,1)==1 && ~iscell(ydat)
            ydat={ydat};
        end
        try
        ydat=ydat(cellfun(@(x) size(x,2)>1,ydat));
        catch
            keyboard
        end
        
        max_val_at_t=NaN(size(ydat,1),size(val,3));
        for lin=1:size(ydat,1)
            if size(ydat{lin},2)==size(val,3)
                max_val_at_t(lin,:)=ydat{lin};
            elseif size(ydat{lin},2)==size(val,3)*9 % each errorbar needs 6 vertices + 3 NaNs
                ttt=reshape(ydat{lin},[9,size(ydat{lin},2)/9]);
                max_val_at_t(lin,:)=max(ttt);
            end
        end
        max_val_at_t=max(max_val_at_t)';
        
        no_trials_not_nan=squeeze(sum(~isnan(val),2));
        %         max_val_at_t=max(val_mean+val_std./sqrt(no_trials_not_nan),[],1)';
        ypos_signmarker=max_val_at_t;
        if no_groups>1
            for group1=1:no_groups
                for group2=group1+1:no_groups
                    
                    ypos_signmarker=ypos_signmarker+double(H(:,group1,group2)==1)*ax_yscale/20;
                    ypos_act=ypos_signmarker;
                    ypos_act(H(:,group1,group2)==0 | isnan(H(:,group1,group2)))=NaN;
                    %                 ypos_signmarker=max_val_at_t+ nansum(nansum(H(:,1:group1,1:group2),2),3)*ax_yscale/20;
                    %                 ypos_signmarker(H(:,group1,group2)==0 | isnan(H(:,group1,group2)))=NaN;
                    sign_handle=plot(1:no_parts,ypos_act,'o','MarkerEdgeColor',colororder(group1,:),'MarkerFaceColor',colororder(group2,:),'MarkerSize',8,'LineWidth',3);
                    set(sign_handle,'UserData',p(:,group1,group2));
                    yl(2)=max(yl(2),max(ypos_act)*1.08);
                end
            end
        end
        set(gca,'YLim',yl);
    end
end