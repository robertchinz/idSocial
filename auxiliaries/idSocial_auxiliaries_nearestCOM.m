function tr_com =idSocial_auxiliaries_nearestCOM(tr)



issctructTr = false;
if iscell(tr) || isstruct(tr)
    issctructTr = true;
    tr_orig = tr;
	tr = idSocial_auxiliaries_formatInputTrajectory(tr);
end
tr = idSocial_auxiliaries_formatInputTrajectory(tr);

no_frames=size(tr,1);
no_focals=size(tr,2);
% if ndims(tr)==3
%     no_dim = size(tr,3);
%     % Add neighbour dimension if necessary
%     tr2=zeros(no_frames,no_fish,no_fish,no_dim);
%     tr2(:,:,:,1:no_dim)=reshape(repmat(tr,[1,no_fish]),no_frames,no_fish,no_fish,no_dim);
%     tr=tr2;
%     clear tr2;
% end    
no_neighbors=size(tr,3);
no_dim = size(tr,4);
%% Calculate distances
foc_to_nb_vec=NaN(no_frames,no_focals,no_neighbors,no_dim);
for ff=1:no_focals
    for nf=1:no_neighbors
        if ff~=nf 
            foc_to_nb_vec(:,ff,nf,:)=squeeze(tr(:,nf,ff,:)-tr(:,ff,ff,:));
        end
    end
end
distances=sqrt(sum(foc_to_nb_vec.^2,4));


[dd nn]=sort(distances(:,:,:),3);
nn(~isfinite(dd))=NaN;
frames_allpresent=sum(sum(isfinite(nn),2),3)==(no_focals*(no_neighbors-1));
%%
% Seperate trajectories from different neighbors by NaNs in
% order to prevent 'strange things' from happening when
% velocities etc. are calculated.
% NOTE: Not necessary, because focal trajectories are not
% affected! 
% nn_change=cat(1,diff(nn,1,1),zeros(1,no_fish,no_fish)) ;
% dd(nn_change~=0)=NaN;
% nn(nn_change~=0)=NaN;
% Only take frames where all individuals are present (if one or more is missing, we cannot be sure if the order is correct!)
nn(~frames_allpresent,:,:)=NaN;
%% Re-order trajectories
% tr_reorder = NaN(size(tr));
% for ff=1:no_focals
%     for nf=1:no_neighbors
%       if nf<no_focals
%         for dim=1:no_dim
%             linidx = sub2ind(size(tr),1:no_frames,squeeze(nn(:,ff,nf))',ones(1,no_frames)*ff,ones(1,no_frames)*dim);
%             linidx_notnan = ~isnan(linidx);
%             tr_reorder(linidx_notnan,nf,ff,dim)=tr(linidx(linidx_notnan));
%         end
%       else
%            tr_reorder(:,nf,ff,:)=tr(:,ff,ff,:);
%             % [nf ff tr(idx,nf,nf,1)]
%       end
%        
%     end
% end
%% Re-order trajectories: Diagonals with orig. trajectories, rest ordered
tr_reorder = NaN(size(tr));
for ff=1:no_focals
    for nf=1:no_neighbors
        if ff ~= nf
            if nf<ff
                for dim=1:no_dim
                    linidx = sub2ind(size(tr),1:no_frames,squeeze(nn(:,ff,nf))',ones(1,no_frames)*ff,ones(1,no_frames)*dim);
                    linidx_notnan = ~isnan(linidx);
                    tr_reorder(linidx_notnan,nf,ff,dim)=tr(linidx(linidx_notnan));
                end
            elseif nf>ff
                for dim=1:no_dim
                    linidx = sub2ind(size(tr),1:no_frames,squeeze(nn(:,ff,nf-1))',ones(1,no_frames)*ff,ones(1,no_frames)*dim);
                    linidx_notnan = ~isnan(linidx);
                    tr_reorder(linidx_notnan,nf,ff,dim)=tr(linidx(linidx_notnan));
                end
            end
        else
            tr_reorder(:,ff,ff,:)=tr(:,ff,ff,:);
            % [nf ff tr(idx,nf,nf,1)]
           
        end
    end
end
%% Calculate coms
tr_com=NaN(no_frames,no_focals,no_neighbors,no_dim);
for ff=1:no_focals
    for nf=1:no_neighbors % Only no_neighbors center of masses
        if ff>nf 
            tr_com(:,nf,ff,:)=squeeze(mean(tr_reorder(:,setxor(ff,1:nf),ff,:),2));
        end
        if ff<nf 
            tr_com(:,nf,ff,:)=squeeze(mean(tr_reorder(:,setxor(ff,1:nf),ff,:),2));
        end
        if ff==nf 
            tr_com(:,nf,ff,:)= tr_reorder(:,nf,ff,:);
        end
    end
end


if issctructTr
    tr_orig.Tr = tr_com;
    tr_com = tr_orig;
end

%%
% dord12 = sqrt((tr_reorder(:,1,1,1)-tr_reorder(:,2,1,1)).^2+(tr_reorder(:,1,1,2)-tr_reorder(:,2,1,2)).^2);
% dord13 = sqrt((tr_reorder(:,1,1,1)-tr_reorder(:,3,1,1)).^2+(tr_reorder(:,1,1,2)-tr_reorder(:,3,1,2)).^2);
% dord14 = sqrt((tr_reorder(:,1,1,1)-tr_reorder(:,4,1,1)).^2+(tr_reorder(:,1,1,2)-tr_reorder(:,4,1,2)).^2);
