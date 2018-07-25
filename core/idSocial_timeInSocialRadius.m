function [timedistribution, mediantime,meantime, time_vals, no_interactions, total_time, total_time_inclNaN, timedistributionRAND, mediantimeRAND,meantimeRAND, time_valsRAND, total_timeRAND, total_time_inclNaNRAND]= ...
    idSocial_timeInSocialRadius(tr,social_radius,edges,bodylength,framerate,random_data)
% Calculates mutual distances between group members

if nargin<2 || isempty(social_radius)
    social_radius=5;
end
if nargin<3 || isempty(edges)
    edges=0:30;
end
if nargin<4 || isempty(bodylength)
    bodylength=1;
end
if nargin<5 || isempty(framerate)
    framerate=1;
end

if nargin<6 || isempty(random_data)
    random_data=false;
end

[tr,~,~,no_frames,no_fish,no_dim] = ...
    idSocial_auxiliaries_formatInputTrajectory(tr);
rand_check = idSocial_auxiliaries_trRandCheck(tr);

foc_to_nb_vec=NaN(no_frames,no_fish,no_fish,2);
for ff=1:no_fish
    for nf=1:no_fish
        if ff~=nf && ~rand_check(ff,nf)
            foc_to_nb_vec(:,ff,nf,:)=squeeze(tr(:,nf,ff,1:2)-tr(:,ff,ff,1:2));
        end
    end
end
distance_focal_neighbour=sqrt(sum(foc_to_nb_vec.^2,4))./bodylength;

distance_focal_neighbour=permute(distance_focal_neighbour,[2,3,1]);
distance_focal_neighbour(logical(repmat(eye(no_fish,no_fish),[1 1 no_frames])))=NaN;
%%
no_bins=size(edges,2)-1;
timedistribution=NaN(no_fish,no_fish,no_bins);
time_vals=cell(no_fish,no_fish);
meantime = NaN(no_fish,no_fish);
mediantime = NaN(no_fish,no_fish);

timedistributionRAND=NaN(no_fish,no_fish,no_bins);
time_valsRAND=cell(no_fish,no_fish);
meantimeRAND = NaN(no_fish,no_fish);
mediantimeRAND = NaN(no_fish,no_fish);
total_time_inclNaN = NaN(no_fish,no_fish);
total_time = NaN(no_fish,no_fish);
no_interactions = NaN(no_fish,no_fish);

total_time_inclNaNRAND = NaN(no_fish,no_fish);
total_timeRAND = NaN(no_fish,no_fish);
for ff=1:no_fish
    for nf=1:no_fish
        if ff~=nf && ~rand_check(ff,nf)
            val=squeeze(distance_focal_neighbour(ff,nf,:));
            %             try
            [onsetIdx, offsetIdx] = idSocial_auxiliaries_findSocialOnOffset(val,social_radius);
            %             catch
            %                 keyboard
            %             end
            
            
            if ~isempty(onsetIdx)
                dt = (offsetIdx-onsetIdx+1)/framerate;
                
                hitemp=histc(dt,edges);
                timedistribution(ff,nf,:) = hitemp(1:end-1);
                time_vals{ff,nf} = dt;
                meantime(ff,nf) = mean(dt);
                mediantime(ff,nf) = median(dt);
                
                no_interactions(ff,nf) = size(onsetIdx,1);
                % Total time closer than social_radius, without pieces
                % containing NaN
                total_time(ff,nf) = sum(dt);
                % Total time, including all frames
                total_time_inclNaN(ff,nf) = sum(val<social_radius)/framerate;
                
                if random_data
                    valShiftRAND = NaN;
                    no_pieces = numel(dt);
                    dtRAND = NaN(no_pieces,1);
                    for pc = 1:no_pieces
                        
                        idx = onsetIdx(pc);
                        trNb = tr(:,nf,ff,:);
                        
                        % Get rid of real cases of interaction
                        
                        trNb(max(1,onsetIdx(pc)-300):min(no_frames,offsetIdx(pc)+300),:,:,:) = NaN;
                        
                        valRAND = sqrt(sum((repmat(tr(idx,ff,ff,1:2),[no_frames,1,1,1])-trNb(:,1,1,1:2)).^2,4))/bodylength;
                        
                        [onsetIdxRAND, offsetIdxRAND] = idSocial_auxiliaries_findSocialOnOffset(valRAND,social_radius);
                        
                        piece_index = [];
                        while isempty(piece_index) && ~isempty(onsetIdxRAND)
                            % Find onset furthest away in time (to get rid of any
                            % actual correlations
                            %                     idxFurthestAway = onsetIdxRAND(max(abs(idx-onsetIdxRAND))==abs(idx-onsetIdxRAND));
                            idxFurthestAway = onsetIdxRAND(randperm(size(onsetIdxRAND,1),1));
                            
                            
                            % Align trajectories of foc and rand-nb
                            trAlign = NaN(no_frames,2,2);
                            trAlign(:,1,:) = tr(:,ff,ff,1:2);
                            
                            if idxFurthestAway>idx
                                shiftIdx = [ (idxFurthestAway-idx+1):no_frames 1:(idxFurthestAway-idx)];
                            elseif idxFurthestAway<idx
                                shiftIdx = [ (idxFurthestAway+(no_frames-idx+1)):no_frames 1:idxFurthestAway+(no_frames-idx)];
                            else
                                shiftIdx = 1:no_frames;
                            end
                            trAlign(:,2,:) = tr(shiftIdx,nf,ff,1:2);
                            
                            % Find time of overlap
                            valShiftRAND = sqrt(sum((trAlign(:,1,1:2)-trAlign(:,2,1:2)).^2,3))/bodylength;
                            
                            [onsetIdxShiftRAND, offsetShiftIdxRAND] = idSocial_auxiliaries_findSocialOnOffset(valShiftRAND,social_radius);
                            
                            piece_index = find(onsetIdxShiftRAND <= idx & offsetShiftIdxRAND >= idx);
                            if ~isempty(piece_index) % According to the construction of the whole thing, there should be at least a 1-frame-piece;
                                % nevertheless if idx lies within the first
                                % overlap, starting at frame 1, it will be
                                % discarded since we do not know when it did
                                % actually start. In this case, piece_index will be
                                % empty.
                                dtRAND(pc) = (offsetShiftIdxRAND(piece_index)-onsetIdxShiftRAND(piece_index)+1)/framerate;
                            else
                                onsetIdxRAND = setxor(onsetIdxRAND,idxFurthestAway);
                            end
                        end
                    end
                    hitempRAND=histc(dtRAND,edges);
                    timedistributionRAND(ff,nf,:) = hitempRAND(1:end-1);
                    time_valsRAND{ff,nf} = dtRAND;
                    meantimeRAND(ff,nf) = mean(dtRAND);
                    mediantimeRAND(ff,nf) = median(dtRAND);
                    total_timeRAND(ff,nf) = sum(dtRAND);
                    if ~isempty(valShiftRAND)
                        total_time_inclNaNRAND(ff,nf) = sum(valShiftRAND<social_radius)/framerate;
                    else
                        total_time_inclNaNRAND(ff,nf) = NaN;
                    end
                    
                end
            end
        end
    end
end

timedistribution=mat2cell(timedistribution,ones(1,no_fish),ones(1,no_fish),ones(1,no_bins));
timedistributionRAND=mat2cell(timedistributionRAND,ones(1,no_fish),ones(1,no_fish),ones(1,no_bins));

end
