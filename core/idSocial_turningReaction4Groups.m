function [turning_mean, turning_median, turning_ratio, turning, funcinfo]= ...
    idSocial_turningReaction4Groups(tr,framerate,bodylength,dist_lim)
% Calculates mutual distances between group members
no_fish=size(tr,2);
no_frames=size(tr,1);
no_dim=size(tr,4);


xax=[1 0 0];    % Setting axis orientation
yax=[0 1 0];    % Setting axis orientation
zax=[0 0 1];    % Setting axis orientation

if nargin<3 || isempty(framerate)
    framerate=2;
end
if nargin<3 || isempty(bodylength)
    bodylength=1;
end

if nargin<4 || isempty(dist_lim)
    dist_lim = [-inf inf];
    
end

invertY = true;
[tr,vel,acc,no_frames,no_fish,no_dim] = ...
    idSocial_auxiliaries_formatInputTrajectory(tr,invertY);
[rand_check, no_fishORIG] = idSocial_auxiliaries_trRandCheck(tr);

[focal_to_neighbour_vector_rotated,vel_rot,focal_a_rotated,turnmatr]= ...
    idSocial_auxiliaries_rotateToFocalSystem(tr, vel, acc);

%%

dim=1;
no_confs = ceil(no_fish/2);

turning_median=NaN(no_fish,no_fish);
turning_mean=NaN(no_fish,no_fish);
turning_ratio=NaN(no_fish,no_fish,no_frames);
turning=cell(no_fish,no_fish,no_confs);
% conf_idces_towards = NaN(no_fish,no_fish,no_frames);
% conf_idces_away = NaN(no_fish,no_fish,no_frames);

no_good_frames = NaN(no_fish,no_fish);
dist = idSocial_distance(tr,bodylength);
xtickLabel = cell(1,no_confs);
xtick = 1:no_confs;

for ff=1:no_fishORIG
    val=focal_a_rotated(:,ff,ff,dim)*framerate^2/bodylength;
    dist_filter = squeeze(dist(ff,:,:))' > dist_lim(1) & squeeze(dist(ff,:,:))' < dist_lim(2);
    acc_right = squeeze(val) > 0 & ~isnan(val);
    acc_left = squeeze(val) < 0 & ~isnan(val);

    for nf=1:no_confs
        if ~rand_check(ff,nf)
            
            conf = [nf-1 no_fishORIG-nf];
            
            xtickLabel{nf} = [num2str(conf(1)) ':' num2str(conf(2))]; 
            
            log_nb_left = squeeze(focal_to_neighbour_vector_rotated(:,:,ff,1))<0;
            log_nb_right = squeeze(focal_to_neighbour_vector_rotated(:,:,ff,1))>0;
            
            allPresent = (sum(log_nb_left,2)  + sum(log_nb_right,2))==no_fishORIG-1; % Seems to be unnecessary: & all(~isnan(squeeze(tr(:,ff,ff,:))),2); 
            
            nb_right = sum(log_nb_right & ...
                ~isnan(squeeze(focal_to_neighbour_vector_rotated(:,:,ff,1))),2) == conf(2) & ...
                sum(dist_filter,2) == no_fish-1;
            nb_left = sum(log_nb_left & ...
                ~isnan(squeeze(focal_to_neighbour_vector_rotated(:,:,ff,1))),2) == conf(2) & ...
                sum(dist_filter,2) == no_fish-1;

            
            ttt = NaN(no_frames,1);
            
            if conf(1)~=conf(2) % If configuration is symmetrical, the next line would render ALL values negative.
                ttt(acc_right & nb_right) = abs(val(acc_right & nb_right));
                ttt(acc_left & nb_left) = abs(val(acc_left & nb_left));
                ttt(acc_right & nb_left | acc_left & nb_right) = -abs(val(acc_right & nb_left | acc_left & nb_right));
            else
                ttt(acc_right) = abs(val(acc_right));
                ttt(acc_left) = -abs(val(acc_left));
            end
            
        
            idxmatr = true(no_frames,no_fishORIG);
            idxmatr(squeeze(isnan(tr(:,:,ff,1)))) = false;
            
       
            acc_right_matr = repmat(acc_right,[1,no_fishORIG]);
            acc_left_matr = repmat(acc_left,[1,no_fishORIG]);
            allPresent_matr = repmat(allPresent,[1,no_fishORIG]);
            
            idces_towards = false(size(log_nb_left));
            idces_towards( ((log_nb_left & acc_left_matr) | ...
                (log_nb_right & acc_right_matr)  )& allPresent_matr) = idxmatr( ((log_nb_left & acc_left_matr) | ...
                (log_nb_right & acc_right_matr)  )& allPresent_matr);
            
            idces_away = false(size(log_nb_left));
            idces_away( ((log_nb_left & acc_right_matr) | ...
                (log_nb_right & acc_left_matr)  )& allPresent_matr) = idxmatr( ((log_nb_left & acc_right_matr) | ...
                (log_nb_right & acc_left_matr)  )& allPresent_matr);
            
%             idces_towards(:,ff) = true(no_frames,1);
%             idces_away(:,ff) = true(no_frames,1);
% 
%             conf_idces_towards(ff,:,:) = idces_towards';
%             conf_idces_away(ff,:,:) = idces_away';
            
            disp(['           Focal ' num2str(ff) ': % for P(' num2str(conf(2)) '|' num2str(conf(1)) ':' num2str(conf(2)) '): ' num2str(sum((acc_right & nb_right | acc_left & nb_left) & ~isnan(ttt))/no_frames)])
            disp(['           Focal ' num2str(ff) ': % for P(' num2str(conf(1)) '|' num2str(conf(1)) ':' num2str(conf(2)) '): ' num2str(sum((acc_right & nb_left | acc_left & nb_right) & ~isnan(ttt))/no_frames)])
            
            ttt = ttt(~isnan(ttt));
            
            turning{ff,ff,no_confs-nf+1}= ttt;
            no_good_frames(ff,no_confs-nf+1) = sum((nb_right | nb_left) & ~isnan(val) & val ~=0 );
          
            turning_mean(ff,no_confs-nf+1)=nanmean(ttt);
            turning_median(ff,no_confs-nf+1)=nanmedian(ttt);
            
            turning_ratio(ff,no_confs-nf+1,1:numel(ttt))=ttt;%idSocial_positive_ratio(ttt,[],0);
                
            
        end
    end
end

fstack=dbstack(3);
funcinfo.Function = mfilename;
funcinfo.callerFunction = fstack(1).name;
funcinfo.no_good_frames = no_good_frames;
funcinfo.no_frames = no_frames;
funcinfo.no_fish = no_fish;
funcinfo.XTick = xtick;
funcinfo.XTickLabel = xtickLabel(end:-1:1);
