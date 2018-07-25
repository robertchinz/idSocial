function [distanceGroupMean, distanceGroupMedian,inbins]=idSocial_groupConfigurationsN4(tr,bodylength,normalization)
% Calculates mutual distances between group members
no_fish=size(tr,2);
no_frames=size(tr,1);
no_dim=size(tr,4);

if nargin<2 || isempty(bodylength)
    bodylength=1;
end
if nargin<3 || isempty(normalization)
    normalization='density';
end


foc_to_nb_vec=NaN(no_frames,no_fish,no_fish,no_dim);
for ff=1:no_fish
    for nf=1:no_fish
        if ff~=nf
            foc_to_nb_vec(:,nf,ff,:)=squeeze(tr(:,nf,ff,:)-tr(:,ff,ff,:));
        end
    end
end
distance_focal_neighbour=sqrt(sum(foc_to_nb_vec.^2,4))./bodylength;
good_frames=sum( ~isnan(distance_focal_neighbour(:,:)),2)==no_fish*(no_fish-1);
if nansum(good_frames)/no_frames<.5
    warning([mfilename ': Only ' num2str( nansum(good_frames)/no_frames*100) '% of frames can be used for calculations.'])
end

distance_focal_neighbour=permute(distance_focal_neighbour,[2,3,1]);
distance_focal_neighbour(logical(repmat(tril(ones(no_fish,no_fish)),[1 1 no_frames ])))=NaN;
distance_focal_neighbour=permute(distance_focal_neighbour,[3,1,2]);

%% 
dist_social_bl = 6;

dmatr = squeeze(distance_focal_neighbour(29,:,:))
dmatr_log = dmatr < dist_social_bl
%% Calculate center of mass for all possible combinations
% 11 combinations (6 + 4 + 1 corresponding to pairs, N=3 and
% N=4
combs=[];%cell(1,1);
for ff=2:no_fish
    combs=vertcat(combs, num2cell(nchoosek(1:no_fish,ff),2));
end

trsq=NaN(no_frames,no_fish, no_dim);
for ff=1:no_fish
    trsq(:,ff,:)=tr(:,ff,ff,:);
end

no_combs = size(combs,1);
com4combs=NaN(no_frames,no_dim,no_combs);
for cb=1:no_combs
    com4combs(:,:,cb)=squeeze(nanmean(trsq(:,combs{cb},:),2));
end

foc_to_com_dist=NaN(no_frames,no_fish,no_combs);
com_to_com_dist=NaN(no_frames,no_combs,no_combs);
for ff=1:no_fish
    for cb=1:no_combs
        
            foc_to_com_dist(:,ff,cb)=sqrt(sum((squeeze(trsq(:,ff,:))-com4combs(:,:,cb)).^2,2))./bodylength;

    end
end
 for cb1=1:no_combs
    for cb2=(1):no_combs
        if cb1 ~=cb2
                    com_to_com_dist(:,cb1,cb2)=sqrt(sum((com4combs(:,:,cb1)-com4combs(:,:,cb2)).^2,2))./bodylength;
        end
    end
end

gr_size=cellfun(@(x) numel(x),combs);
member_filter =  cellfun(@(x) ismember(1:no_fish,x),combs,'Uniformoutput', false);
member_filter =  vertcat(member_filter{:})';
%%
testvec=squeeze(foc_to_com_dist(1000,:,:))<dist_social_bl

testvec(~member_filter)=false; % Only those indiv. from which the com has been calculated should be in the group

% Check for which group indiv. are close enough to their
% com:
com_to_com_dist(:,ff,cb)
good_groups = sum(testvec,1)==gr_size';


% Check if the coms of the found groups are far apart (there can be more than one group at the same time!):
comtestvec = squeeze(com_to_com_dist(1000,:,:))>2*dist_social_bl;
% comtestvec(logical(eye(size(comtestvec))))=true;
very_good_groups = logical(double(good_groups)*double(comtestvec));

% gg = find(good_groups);
% [com1 com2]=find(comtestvec);
% intersect(gg,[com1 com2],'rows');
sort_of_group = gr_size(very_good_groups);
%%
% 
% cmb_act=11;
% figure; 
% plot(squeeze(tr(29,1,1,1)),squeeze(tr(29,1,1,2)),'square')
% % axis equal
% hold on
% plot(squeeze(tr(29,2,2,1)),squeeze(tr(29,2,2,2)),'square')
% plot(squeeze(tr(29,3,3,1)),squeeze(tr(29,3,3,2)),'square')
% plot(squeeze(tr(29,4,4,1)),squeeze(tr(29,4,4,2)),'square')
% plot(com4combs(29,1,cmb_act),com4combs(29,2,cmb_act),'r*')

%%

%%
% if   ~any(dmatr_log(:)) % All indiv. 'asocial'
%     
% elseif  % 2 + 2
% elseif % 1 + 3
% elseif % 4 + 0
% end
%%
keyboard
distanceGroupMean=NaN;
distanceGroupMedian=NaN;
inbins=cell(no_fish,no_fish,no_fish*(no_fish-1)/2);

[sortedGrDist, sortedGrDistIdx] =sort(distance_focal_neighbour(:,:),2);
sortedGrDist(~good_frames,:)=NaN;
sortedGrDist=sortedGrDist(:,1:no_fish*(no_fish-1)/2);

distanceGroupMean = nanmean(sortedGrDist,1);

inbins(1,1,:)=num2cell(sortedGrDist(good_frames,:),1);
% %%
% coordinates=sortedGrDist(:,[1 2]);
% edges={0:1:10;0:1:10};
% [heatmap,...
%                     ~, ~, maps_temp,...
%                     stdmaps_temp,...
%                     indices_in_bin]=...
%                     idSocial_hist3indexmap_mosaic(coordinates,'map',[],'edges',edges,'spacing',1);
%                  figure; imagesc(edges{1},edges{2},bsxfun(@rdivide,heatmap,sum(heatmap,1))); set(gca,'YDir','normal')