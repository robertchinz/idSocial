function varargout = ...
    idSocial_timeInSocialRadiusRank(tr,social_radius,no_points,bodylength,framerate,random_data)
% Calculates mutual distances between group members

[tr,~,~,no_frames,no_fish,no_dim] = ...
    idSocial_auxiliaries_formatInputTrajectory(tr);
rand_check = idSocial_auxiliaries_trRandCheck(tr);

if nargin<2 || isempty(social_radius)
    social_radius=5;
end
if nargin<3 || isempty(no_points)
    no_points=100;
end
if nargin<4 || isempty(bodylength)
    bodylength=1;
end
if nargin<5 || isempty(framerate)
    framerate=1;
end

if nargin<6 || isempty(random_data) || strcmpi(random_data,'false')
    random_data=[];
    no_neighborsRAND = no_fish;
elseif ~isempty(random_data)
    if isnumeric(random_data)
        no_neighborsRAND = random_data;
    elseif islogical(random_data)
        no_neighborsRAND = no_fish;
    end
end

no_neighborsRAND_Orig = no_neighborsRAND;
if no_fish >2
    no_neighborsRAND = no_fish;
end


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
% no_bins=size(edges,2)-1;
no_bins = no_points;
timedistribution=NaN(no_fish,no_fish,no_bins);
time_vals=cell(no_fish,no_fish);
meantime = NaN(no_fish,no_fish);
mediantime = NaN(no_fish,no_fish);
    
timedistributionRAND=NaN(no_fish,no_neighborsRAND_Orig,no_bins);
time_valsRAND=cell(no_fish,no_neighborsRAND_Orig);
meantimeRAND = NaN(no_fish,no_neighborsRAND_Orig);
mediantimeRAND = NaN(no_fish,no_neighborsRAND_Orig);
total_time_inclNaN = NaN(no_fish,no_fish);
total_time = NaN(no_fish,no_fish);
no_interactions = NaN(no_fish,no_fish);
no_interactionsAll = NaN(no_fish,no_fish);

total_time_inclNaNRAND = NaN(no_fish,no_neighborsRAND_Orig);
total_timeRAND = NaN(no_fish,no_neighborsRAND_Orig);

onsetIdxCell = cell(no_fish,no_fish);
for ff=1:no_fish
    onsetIdx = cell(1,no_fish);
    offsetIdx = cell(1,no_fish);
    for nf=1:no_fish
        if ff~=nf && ~rand_check(ff,nf)
            val=squeeze(distance_focal_neighbour(ff,nf,:));
            %             try
            [onsetIdx{nf}, offsetIdx{nf}] = idSocial_auxiliaries_findSocialOnOffset(val,social_radius);
            %             catch
            %                 keyboard
            %             end
            onsetIdxCell{ff,nf} = onsetIdx{nf};
            
            no_interactionsAll(ff,nf) = size(onsetIdx{nf},1);
            if ~isempty(onsetIdx{nf}) && no_fish<=2
                dt = (offsetIdx{nf}-onsetIdx{nf}+1)/framerate;
                
                %                 hitemp=histc(dt,edges);
                hitemp = sort(dt(~isnan(dt)),'descend');
                if ~isempty(hitemp)
                    timedistribution(ff,nf,1:min(numel(hitemp),no_bins)) = hitemp(1:min(numel(hitemp),no_bins));%(1:end-1);
                end
                time_vals{ff,nf} = hitemp;
                meantime(ff,nf) = mean(dt);
                mediantime(ff,nf) = median(dt);
                
                no_interactions(ff,nf) = size(onsetIdx{nf},1);
                
                
                % Total time closer than social_radius, without pieces
                % containing NaN
                total_time(ff,nf) = sum(dt);
                % Total time, including all frames
                total_time_inclNaN(ff,nf) = sum(val<social_radius)/framerate;
            end
        end
    end
    
    if no_fish>2  % For N>2: Combine overlapping neighbors. Distiction between individual neighbors is lost!
        
        allNbsIdces = [vertcat(onsetIdx{:}) vertcat(offsetIdx{:})];
        if ~isempty(allNbsIdces)
            allNbsIdcesSort = sortrows(allNbsIdces,1);
            no_idces = size(allNbsIdcesSort,1);
            count=1;
            while count<no_idces &&  size(allNbsIdcesSort,1)>=count
                %
                includedIdx = allNbsIdcesSort(:,2)<allNbsIdcesSort(count,2) & allNbsIdcesSort(:,1)>allNbsIdcesSort(count,1) ;
                if any(includedIdx)
                    allNbsIdcesSort(includedIdx,:)=[];
                else
                    count=count+1;
                end
                
            end
            count=1;
            no_idces = size(allNbsIdcesSort,1);
            while count<no_idces &&  size(allNbsIdcesSort,1)>=count
                overlapIdx = allNbsIdcesSort(:,1)>allNbsIdcesSort(count,1) & allNbsIdcesSort(:,2)>allNbsIdcesSort(count,2) & ...
                    allNbsIdcesSort(:,1)<allNbsIdcesSort(count,2);
                if any(overlapIdx)
                    allNbsIdcesSort(count,1) = allNbsIdcesSort(count,1);
                    allNbsIdcesSort(count,2) = max(allNbsIdcesSort(overlapIdx,2));
                    allNbsIdcesSort(overlapIdx,:) = [];
                else
                    count=count+1;
                end
            end
            
            
            if ~isempty(allNbsIdcesSort)
                
                nf=1;
                dt = (allNbsIdcesSort(:,2)-allNbsIdcesSort(:,1)+1)/framerate;
                
                %                 hitemp=histc(dt,edges);
                hitemp = sort(dt(~isnan(dt)),'descend');
                if ~isempty(hitemp)
                    timedistribution(ff,nf,1:min(numel(hitemp),no_bins)) = hitemp(1:min(numel(hitemp),no_bins));%(1:end-1);
                end
                time_vals{ff,nf} = hitemp;
                meantime(ff,nf) = mean(dt);
                mediantime(ff,nf) = median(dt);
                
                no_interactions(ff,nf) = size(allNbsIdcesSort,1);
                % Total time closer than social_radius, without pieces
                % containing NaN
                total_time(ff,nf) = sum(dt);
                % Total time, including all frames
                total_time_inclNaN(ff,nf) = NaN;%sum(val<social_radius)/framerate;
            end
        end
    end
    
    if ~isempty(random_data)
        count2=1;
        onsetIdxShiftRAND2 = cell(1,no_neighborsRAND*max(no_interactionsAll(:)));
        possible_nbs = setxor(ff,1:no_fish);
        possible_indices = find(~any(isnan(tr(:,ff,ff,:)),4)); % Exclude focal indices where coordinates are NaN.
        for nfRAND=1:no_neighborsRAND
            if ff~=nfRAND && ~rand_check(ff,nf)
                nf = possible_nbs(mod(nfRAND-1,numel(possible_nbs))+1);%   mod(nfRAND-1,no_fish)+1;
                
                no_pieces = no_interactionsAll(ff,nf);
                if ~isnan(no_pieces) && no_pieces>0
                    dtRAND = NaN(no_pieces,1);
                    randIdx = possible_indices(randi(numel(possible_indices),1,no_pieces));
                    for pc = 1:no_pieces
                        
                        idx = randIdx(pc); % Instead of using the focal's original onset indices,
                        % I also randomize these because if
                        % no_neighborsRAND > no_fish, the same interactions will occur for various
                        % random neighbors because the number of
                        % encounters for (idx-of-focal-at-real-onset, nb-at-position-of-focal) is limited
                        %                     idx = onsetIdxCell{ff,nf}(pc);
                        trNb = tr(:,nf,ff,:);
                        
                        % Get rid of real cases of interaction
                        
                        trNb(max(1,idx-300):min(no_frames,idx+300),:,:,:) = NaN;
                        
                        valRAND = sqrt(sum((repmat(tr(idx,ff,ff,1:2),[no_frames,1,1,1])-trNb(:,1,1,1:2)).^2,4))/bodylength;
                        
                        [onsetIdxRAND, offsetIdxRAND] = idSocial_auxiliaries_findSocialOnOffset(valRAND,social_radius);
                        
                        piece_index = [];
                        valShiftRAND = [];
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
                                onsetIdxShiftRAND2{count2} = [onsetIdxShiftRAND(piece_index) offsetShiftIdxRAND(piece_index)];
                                count2 = count2+1;
                            else
                                onsetIdxRAND = setxor(onsetIdxRAND,idxFurthestAway);
                            end
                        end
                    end
                    %                     hitempRAND=histc(dtRAND,edges);
                    hitempRAND = sort(dtRAND(~isnan(dtRAND)),'descend');
                    %                     hitempRAND = hitempRAND(1:no_points);
                    if ~isempty(hitempRAND)
                        timedistributionRAND(ff,nfRAND,1:min(numel(hitempRAND),no_bins)) = hitempRAND(1:min(numel(hitempRAND),no_bins));%(1:end-1);
                    end
                    if no_fish<=2
                        time_valsRAND{ff,nfRAND} = hitempRAND;
                        meantimeRAND(ff,nfRAND) = mean(dtRAND);
                        mediantimeRAND(ff,nfRAND) = median(dtRAND);
                        total_timeRAND(ff,nfRAND) = sum(dtRAND);
                        if ~isempty(valShiftRAND)
                            total_time_inclNaNRAND(ff,nfRAND) = sum(valShiftRAND<social_radius)/framerate;
                        end
                    end
                end
            end
        end
        
        if no_fish>2  % For N>2: Combine overlapping neighbors. Distiction between individual neighbors is lost!
            
            allNbsIdces = vertcat(onsetIdxShiftRAND2{:});
            if ~isempty(allNbsIdces)
                allNbsIdcesSort = sortrows(allNbsIdces,1);
                no_idces = size(allNbsIdcesSort,1);
                count=1;
                while count<no_idces &&  size(allNbsIdcesSort,1)>=count
                    %
                    includedIdx = allNbsIdcesSort(:,2)<allNbsIdcesSort(count,2) & allNbsIdcesSort(:,1)>allNbsIdcesSort(count,1) ;
                    if any(includedIdx)
                        allNbsIdcesSort(includedIdx,:)=[];
                    else
                        count=count+1;
                    end
                    
                end
                count=1;
                no_idces = size(allNbsIdcesSort,1);
                while count<no_idces &&  size(allNbsIdcesSort,1)>=count
                    overlapIdx = allNbsIdcesSort(:,1)>allNbsIdcesSort(count,1) & allNbsIdcesSort(:,2)>allNbsIdcesSort(count,2) & ...
                        allNbsIdcesSort(:,1)<allNbsIdcesSort(count,2);
                    if any(overlapIdx)
                        allNbsIdcesSort(count,1) = allNbsIdcesSort(count,1);
                        allNbsIdcesSort(count,2) = max(allNbsIdcesSort(overlapIdx,2));
                        allNbsIdcesSort(overlapIdx,:) = [];
                    else
                        count=count+1;
                    end
                end
                
                
                if ~isempty(allNbsIdcesSort)
                    
                    nf=1;
                    dt = (allNbsIdcesSort(:,2)-allNbsIdcesSort(:,1)+1)/framerate;
                    
                    %                 hitemp=histc(dt,edges);
                    hitempRAND = sort(dt(~isnan(dt)),'descend');
                    if ~isempty(hitemp)
                        timedistribution(ff,nf,1:min(numel(hitemp),no_bins)) = hitemp(1:min(numel(hitemp),no_bins));%(1:end-1);
                    end
                    time_valsRAND{ff,nf} = hitempRAND;
                    meantimeRAND(ff,nf) = mean(dt);
                    mediantimeRAND(ff,nf) = median(dt);
                    
                    % Total time closer than social_radius, without pieces
                    % containing NaN
                    total_timeRAND(ff,nf) = sum(dt);
                    % Total time, including all frames
                    total_time_inclNaNRAND(ff,nf) = NaN;%sum(valShiftRAND<social_radius)/framerate;
                end
            end
        end
        
        
        
        
    end
end
timedistribution=mat2cell(timedistribution,ones(1,no_fish),ones(1,no_fish),ones(1,no_bins));
varargout{1} = timedistribution;
varargout{2} = mediantime;
varargout{3} = meantime;
varargout{4} = time_vals;
varargout{5} = no_interactions;
varargout{6} = total_time;
varargout{7} = total_time_inclNaN;

% if ~isempty(random_data)
timedistributionRAND=mat2cell(timedistributionRAND,ones(1,no_fish),ones(1,no_neighborsRAND_Orig),ones(1,no_bins));
varargout{8} = timedistributionRAND;
varargout{9} = mediantimeRAND;
varargout{10} = meantimeRAND;
varargout{11} = time_valsRAND;
varargout{12} = NaN(no_fish,no_neighborsRAND_Orig);
varargout{13} = total_timeRAND;
varargout{14} = total_time_inclNaNRAND;
% end


end
