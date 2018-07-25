function [turningReact, turningReactInBins, turningReactCell,turningReactOneBin]=idSocial_attractionReaction(tr,spacing,edges,framerate,bodylength,normalize_by_individual_max,rotate_z,method)

no_fish=size(tr,2);
no_frames=size(tr,1);
no_dim=size(tr,4);

xax=[1 0 0];    % Setting axis orientation
yax=[0 1 0];    % Setting axis orientation
zax=[0 0 1];    % Setting axis orientation

if nargin<4 || isempty(framerate)
    framerate=1;
end
if nargin<5 || isempty(bodylength)
    bodylength=1;
end

if nargin<6 || isempty(normalize_by_individual_max)
    normalize_by_individual_max=false;
end
if nargin<7 || isempty(rotate_z)
    rotate_z=false;
end
if nargin<8 || isempty(method)
    method='hist';
end
if iscell(edges)
    edges=edges{1};
end

invertY = true;
[tr,vel,acc,no_frames,no_fish,no_dim] = ...
    idSocial_auxiliaries_formatInputTrajectory(tr,invertY);
rand_check = idSocial_auxiliaries_trRandCheck(tr);

foc_to_nb_vec=NaN(no_frames,no_fish,no_fish,3);
vel_magn=sqrt(sum(vel.^2,4));
vel_norm=bsxfun(@rdivide,vel,vel_magn);

for ff=1:no_fish
    for nf=1:no_fish
        if ff~=nf && ~rand_check(ff,nf)
            foc_to_nb_vec(:,nf,ff,:)=squeeze(tr(:,nf,ff,:)-tr(:,ff,ff,:));
        end
    end
end


foc_to_nb_distance = sqrt(sum(foc_to_nb_vec.^2,4));
foc_to_nb_dir=foc_to_nb_vec./repmat(foc_to_nb_distance,[1,1,1,3]);

acc_dir=acc./repmat(sqrt(sum(acc.^2,4)),[1,1,1,3]);

[focal_to_neighbour_vector_rotated,focal_v_rotated,focal_a_rotated] = ...
    idSocial_auxiliaries_rotateToFocalSystem(tr,vel,acc);


%%
dim=1;
no_bins=size(edges,2)-1;
no_fish=size(focal_to_neighbour_vector_rotated,2);

turningReactOneBin=NaN(no_fish,no_fish,no_frames);
turningReact=NaN(no_fish,no_fish,no_bins);
turningReactInBins=cell(no_fish,no_fish,no_bins);
for f1=1:no_fish
    for f2=1:no_fish
        if f1~=f2
            
            distance_vector=foc_to_nb_distance(:,f2,f1)/bodylength;
%             map_vector=focal_a_rotated(:,f1,f1,dim)*framerate^2/bodylength;
            map_vector=(sum(squeeze(acc_dir(:,f1,f1,:).*foc_to_nb_dir(:,f2,f1,:)),2));
            %             figure; hist(map_vector(:),200)
            %             map_vector(abs(map_vector)<50)=NaN;
            %             disp('del this!')
          
            if normalize_by_individual_max
                map_vector=map_vector/max(abs(map_vector));
            end
            
            
            %%%%
            %             spacing=1;
            if strcmp(method,'gauss')
                ctrs=edges(1:end-1)+(edges(2)-edges(1))/2;
                prob=exp(-.5*((repmat(distance_vector',[1 numel(ctrs)])-repmat(ctrs,[no_frames 1]))/spacing).^2);
                A = prob.*repmat(map_vector,[1 numel(ctrs)]);
                
                turningReact=nanmean(A,1);
                
                turningReactInBins(f1,f2,:)=arrayfun(@(k) A(~isnan(A(:,k)),k),1:numel(ctrs),'UniformOutput',false);
                
                %             figure; plot(nansum(A,1))
                %%%%
            else
                try
                    [~, bin]=histc(distance_vector,edges);
                catch
                    keyboard
                end
                
                good_idx=bin~=0 & bin<size(edges,2);
                
                lag=spacing;
                
                A=accumarray(bin(good_idx),map_vector(good_idx),[size(edges,2)-1 1],@nanmean);
                
                
                A=smooth(A,'moving',lag);
                
                turningReact(f1,f2,ceil(lag/2):end-floor(lag/2))=A(ceil(lag/2):end-floor(lag/2));
                
                B=accumarray(bin(good_idx),map_vector(good_idx),[size(edges,2)-1 1] ,@(x) {x},{NaN});
                
                
                turningReactInBins(f1,f2,ceil(lag/2):end-floor(lag/2))=B(ceil(lag/2):end-floor(lag/2));
            end
            %             keyboard
            %%%%%
            ttt = map_vector.*distance_vector;
            ttt = ttt(~isnan(ttt));
            turningReactOneBin(f1,f2,1:numel(ttt))= ttt>0;
        end
    end
end

turningReactCell=num2cell(turningReact);