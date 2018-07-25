function [H, p, CI, stat, sampstat]= ...
    idSocial_auxiliaries_statisticalSignificance(val,statistical_test_type,options) %,statistical_test_significance,bootstrap_repetitions,display_mode,bootstrap_statfunc,tail,compare_rand_only,mode_edges,mode_bandwidth,mode_method,sameSampleSize)
% Fomat val: (groups,trials,parts)
test_list = {'ranksum', ...
    'polaviejaStat', ...
    'polaviejaStatRandControls', ...
    'ttest2'};
if nargin < 2
    H = test_list;
    return;
end

no_groups=size(val,1);
no_parts = size(val,3);
no_parts2 = size(val,4);



if nargin<2 || isempty(statistical_test_type)
    statistical_test_type='ranksum';
end


switch statistical_test_type
    case 'ranksum'
        def_options.alpha = 0.05;
        def_options.method = {'default' 'exact' 'approximate'};
        def_options.tail = {'both', 'right', 'left'};
    case {'polaviejaStat','polaviejaStatRandControls'}
        def_options.alpha = 0.05;
        def_options.method = {'mean','median','var','positive_ratio'};
        def_options.tail = {'both', 'right', 'left'};
        def_options.repetitions = 1000;
        def_options.sameSampleSize = true;
        
        
    case 'ttest2'
        def_options.m = 0;
        def_options.alpha = 0.05;
        def_options.Dim = 'default';
        def_options.tail = {'both', 'right', 'left'};
        def_options.Vartype = {'equal','unequal'};
end


if nargin == 2
    H = def_options;
    return;
end

def_options.one_way_only = false;

[~, options_new]=idSocial_readparams([],options,def_options,[]);
if ~isfield(def_options,'repetitions')
    reps=1;
else
    reps = options_new.repetitions;
end
stat=NaN(no_groups,no_groups,no_parts,no_parts2,reps);
%% Calculate statistical significance between groups


CI=[];
H=NaN(no_groups,no_groups,no_parts,no_parts2);
p=NaN(no_groups,no_groups,no_parts,no_parts2);
sampstat=NaN(no_groups,no_groups,no_parts,no_parts2);
if options_new.one_way_only
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


for act_idx=1:size(group_index,1)
    %     for group2=group_index
    group=group_index(act_idx,1);
    group2=group_index(act_idx,2);
    for step2 = 1:no_parts2
        if group~=group2
            switch statistical_test_type
                case 'ttest2'
                    if strcmpi(options_new.Dim,'default')
                        [H(group,group2,:,step2), p(group,group2,:,step2)]=ttest2(squeeze(val(group,:,:,step2)),squeeze(val(group2,:,:,step2)),'Alpha',options_new.alpha,'Tail',options_new.tail,'Vartype',options_new.Vartype);
                    else
                        [H(group,group2,:,step2), p(group,group2,:,step2)]=ttest2(squeeze(val(group,:,:,step2)),squeeze(val(group2,:,:,step2)),'Alpha',options_new.alpha,'Tail',options_new.tail,'Vartype',options_new.Vartype,'Dim',options_new.Dim);
                    end
                case 'bootstrap'
                    for step=1:no_parts
                        [H(group,group2,step,step2), p(group,group2,step,step2), CI, stat(group,group2,step,step2,:), sampstat(group,group2,step,step2)]= ...
                            idSocial_bootstrap(squeeze(val(group,:,step,step2))',squeeze(val(group2,:,step,step2))',reps,options_new.alpha,options_new.method, options_new.tail);
                    end
                case 'polaviejaStat'
                    
                    for step=1:no_parts
                        
                        try
                            
                            [H(group,group2,step,step2), p(group,group2,step,step2), CI, stat(group,group2,step,step2,:), sampstat(group,group2,step,step2)]= ...
                                idSocial_polaviejaStat(squeeze(val(group,:,step,step2)),squeeze(val(group2,:,step,step2)),reps,options_new.alpha,...
                                options_new.method, options_new.tail, [], [], [],options_new.sameSampleSize);
                            
                        catch
                            keyboard
                        end
                    end
                case 'polaviejaStatRandControls'
                    
                    for step=1:no_parts
                        
                        
                        [H(group,group2,step,step2), p(group,group2,step,step2), CI, stat(group,group2,step,step2,:), sampstat(group,group2,step,step2)]= ...
                            idSocial_polaviejaStat4RandomizedControls(squeeze(val(group,:,step,step2))',squeeze(val(group2,:,step,step2))',reps,options_new.alpha,...
                            options_new.method, options_new.tail, [], [], []);
                        
                    end
                case 'kstest2'
                    for step=1:no_parts
                        [p(group,group2,step,step2), H(group,group2,step,step2)]=kstest2(squeeze(val(group,~isnan(val(group,:,step,step2)),step))',squeeze(val(group2,~isnan(val(group2,:,step,step2)),step)),options_new.alpha);
                    end
                otherwise
                    statistical_test_type='ranksum';
                    vr = version('-release');
                    vr_yy = str2double(vr(1:4));
                    vr_aa = vr(end);
                    if (vr_yy == 2012 && vr_aa > 'a') || vr_yy > 2012
                        for step=1:no_parts
                            if ~(all(isnan(val(group,:,step,step2)))) && ~(all(isnan(val(group2,:,step,step2))))
                                if strcmp(options_new.method,'default')
                                    [p(group,group2,step,step2), H(group,group2,step,step2)]=ranksum(squeeze(val(group,~isnan(val(group,:,step,step2)),step))',squeeze(val(group2,~isnan(val(group2,:,step,step2)),step)),'alpha',options_new.alpha,'tail',options_new.tail);
                                else
                                    [p(group,group2,step,step2), H(group,group2,step,step2)]=ranksum(squeeze(val(group,~isnan(val(group,:,step,step2)),step))',squeeze(val(group2,~isnan(val(group2,:,step,step2)),step)),'alpha',options_new.alpha,'tail',options_new.tail, ...
                                        'method',options_new.method);
                                end
                            end
                        end
                    else
                        if ~strcmp(options_new.tail,'both')
                            warning([mfilename ': No tail options_new available for ranksum in Matlab versions < 2012b. Using both sided test!'])
                            disp([mfilename ': No tail options_new available for ranksum in Matlab versions < 2012b. Using both sided test!'])
                        end
                        for step=1:no_parts
                            if ~(all(isnan(val(group,:,step,step2)))) && ~(all(isnan(val(group2,:,step,step2))))
                                if strcmp(options_new.method,'default')
                                    [p(group,group2,step,step2), H(group,group2,step,step2)]=ranksum(squeeze(val(group,~isnan(val(group,:,step,step2)),step))',squeeze(val(group2,~isnan(val(group2,:,step,step2)),step)),'alpha',options_new.alpha);
                                else
                                    [p(group,group2,step,step2), H(group,group2,step,step2)]=ranksum(squeeze(val(group,~isnan(val(group,:,step,step2)),step))',squeeze(val(group2,~isnan(val(group2,:,step,step2)),step)),'alpha',options_new.alpha,...
                                        'method',options_new.method);
                                end
                            end
                        end
                    end
                    
            end
            
            %     end
        end
    end
end

% if nargin>3 && ~isempty(ax) && ishandle(ax)
%
%     val_mean=squeeze(nanmean(val,2));
%     val_std=squeeze(nanstd(val,[],2));
%
%     yl=get(ax,'YLim');
%
%
%     ax_yscale=diff(get(ax,'YLim'));
%
%
%     if strcmp(display_mode,'hist')
%
%         lineObj=findobj(get(ax,'Children'),'Type','line','-or','Type','stair');
%         ydat=get(lineObj,'YData');
%         if size(ydat,1)==1 && ~iscell(ydat)
%             ydat={ydat};
%         end
%         ydat=ydat(cellfun(@(x) size(x,2)>1,ydat));
%         pvaluefontsize=10;
%         linewidth=2;
%         for part=1:size(H,1)
%             for group1=1:no_groups
%                 for group2=group1+1:no_groups
%                     max_val_at_t=max(ydat{group1},[],2)';
%                     yrange=yl(2)*.9-max_val_at_t;
% %                     ypos_signmarker=max_val_at_t;
%
%
%                     ypos_signmarker=max_val_at_t+nansum(nansum(nansum(squeeze(H(1:part,1:group1,1:group2)),1),2))/(nansum(H(:))/2)*yrange;
%                     ypos_signmarker(H(part,group1,group2)==0 | isnan(H(part,group1,group2)))=NaN;
%                     text((nanmean(val(group1,:,part))+nanmean(val(group2,:,part)))/2,ypos_signmarker,['p=' num2str(p(part,group1,group2),'%.2g')],...
%                         'Color',colororder(group2,:),'HorizontalAlignment','center','VerticalAlignment','bottom','fontsize',pvaluefontsize,...
%                         'BackgroundColor','w','Margin',1)
%
%                     plot([nanmean(val(group1,:,part)) nanmean(val(group2,:,part))],[ypos_signmarker ypos_signmarker],'-','Color',colororder(group1,:),'MarkerFaceColor',colororder(group2,:),'MarkerSize',8,'LineWidth',linewidth)
%                     plot([nanmean(val(group1,:,part)) nanmean(val(group1,:,part))],yl,'-','Color',colororder(group1,:),'MarkerFaceColor',colororder(group1,:),'MarkerSize',8,'LineWidth',linewidth)
%                     plot([nanmean(val(group2,:,part)) nanmean(val(group2,:,part))],yl,'-','Color',colororder(group2,:),'MarkerFaceColor',colororder(group2,:),'MarkerSize',8,'LineWidth',linewidth)
%
%                 end
%             end
%         end
%
%
%     else
%
%         lineObj=findobj(get(ax,'Children'),'Type','line');
%         if isempty(lineObj)
%          lineObj=findobj(get(ax,'Children'),'Type','errorbar');
%         end
%         if isfield(get(lineObj(1)),'UData')
%              ydat=cellfun(@(x,y) x+y,get(lineObj,'YData'),get(lineObj,'UData'),'UniformOutput',false);
%         else
%             ydat=get(lineObj,'YData');
%         end
%         if size(ydat,1)==1 && ~iscell(ydat)
%             ydat={ydat};
%         end
%         try
%         ydat=ydat(cellfun(@(x) size(x,2)>1,ydat));
%         catch
%             keyboard
%         end
%
%         max_val_at_t=NaN(size(ydat,1),size(val,3));
%         for lin=1:size(ydat,1)
%             if size(ydat{lin},2)==size(val,3)
%                 max_val_at_t(lin,:)=ydat{lin};
%             elseif size(ydat{lin},2)==size(val,3)*9 % each errorbar needs 6 vertices + 3 NaNs
%                 ttt=reshape(ydat{lin},[9,size(ydat{lin},2)/9]);
%                 max_val_at_t(lin,:)=max(ttt);
%             end
%         end
%         max_val_at_t=max(max_val_at_t)';
%
%         no_trials_not_nan=squeeze(sum(~isnan(val),2));
%         %         max_val_at_t=max(val_mean+val_std./sqrt(no_trials_not_nan),[],1)';
%         ypos_signmarker=max_val_at_t;
%         if no_groups>1
%             for group1=1:no_groups
%                 for group2=group1+1:no_groups
%
%                     ypos_signmarker=ypos_signmarker+double(H(:,group1,group2)==1)*ax_yscale/20;
%                     ypos_act=ypos_signmarker;
%                     ypos_act(H(:,group1,group2)==0 | isnan(H(:,group1,group2)))=NaN;
%                     %                 ypos_signmarker=max_val_at_t+ nansum(nansum(H(:,1:group1,1:group2),2),3)*ax_yscale/20;
%                     %                 ypos_signmarker(H(:,group1,group2)==0 | isnan(H(:,group1,group2)))=NaN;
%                     sign_handle=plot(1:no_parts,ypos_act,'o','MarkerEdgeColor',colororder(group1,:),'MarkerFaceColor',colororder(group2,:),'MarkerSize',8,'LineWidth',3);
%                     set(sign_handle,'UserData',p(:,group1,group2));
%                     yl(2)=max(yl(2),max(ypos_act)*1.08);
%                 end
%             end
%         end
%         set(gca,'YLim',yl);
%     end
% end