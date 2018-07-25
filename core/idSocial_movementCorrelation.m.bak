function [maxcorr, delay_rel,delay_rel_abs,means, orden_out, correlmedia]=idSocial_movementCorrelation(trayectorias,max_delay_frames,max_distance_bl,delay_correlation_lim,framerate,blpxl)
%   Pairwise correlation of movement direction.
%
%   maxcorr = ...
%           IDSOCIAL_MOVEMENTCORRELATION(TR) 
%   calculates the average maximum correlation coefficient.   
%
%   See also idSocial_plotGraph, idSocial_delay2grafo 

%   2014 Robert C. Hinz, Gonzalo G. de Polavieja, Consejo Superior de Investigaciones Científicas

if nargin < 2 || isempty(max_delay_frames)
    max_delay_frames = 60;
end
if nargin < 3 || isempty(max_distance_bl)
    max_distance_bl=inf;
end
if nargin < 4 || isempty(delay_correlation_lim)
    delay_correlation_lim=.5;
end
if nargin < 5 || isempty(framerate)
    framerate=1;
end

max_distance_in_pxl = max_distance_bl * blpxl;

if ndims(trayectorias) == 4
    trayectorias=squeeze(nanmean(trayectorias,3));
end
no_fish=size(trayectorias,2);
no_frames=size(trayectorias,1);


correl_delay=NaN(no_frames-1,max_delay_frames+1);
% delayed_dist=NaN(no_frames-1,no_fish,no_fish,max_delay_frames+1);
correlmedia=NaN(no_fish,no_fish,2*max_delay_frames+1);
correlmedia_noPts = NaN(no_fish,no_fish);

% Ahora correlaciones
delay=-max_delay_frames:1:max_delay_frames;
% maxcorr(1:no_fish,1:no_fish)=0;
% delay_rel(1:no_fish,1:no_fish)=0;
maxcorr= zeros(no_fish,no_fish);
delay_rel = zeros(no_fish,no_fish);


for c1_peces=1:no_fish
    for c2_peces=c1_peces+1:no_fish
        dist=sqrt(sum((trayectorias(:,c1_peces,:)-trayectorias(:,c2_peces,:)).^2,3));
        frlist = dist<max_distance_in_pxl;
%         if any(frlist)
            for c_delay=1:length(delay)
                [correl_act, ~, delayed_dist_act]=trayectorias2correlacion(trayectorias(:,[c1_peces c2_peces],:),delay(c_delay));
                correl_delay(:,c_delay)=correl_act(:,1,2);
                delayed_dist=delayed_dist_act(:,1,2);
            end
            frlist = (delayed_dist<max_distance_in_pxl);
            correlmedia(c1_peces,c2_peces,:)=nanmean(correl_delay(frlist(1:end-1),:));
            correlmedia_noPts(c1_peces,c2_peces) = sum(frlist);
            idx = find(frlist(1:end-1));
            figure;
            plot(delay,squeeze(correl_delay(idx(1:60),:)));
            hold on
            plot(delay,squeeze(correlmedia(c1_peces,c2_peces,:)),'LineWidth',2);
            hold on;
            
            correlmedia_noPts = sum(isfinite(correl_delay(frlist(1:end-1),:)),2);
            [maxcorr(c1_peces,c2_peces),ind]=max(correlmedia(c1_peces,c2_peces,:));
            
            if ~(ind == 1 || ind == size(correlmedia,3))
                delay_rel(c1_peces,c2_peces)=delay(ind)/framerate;
            else
                delay_rel(c1_peces,c2_peces) = NaN;
                maxcorr(c1_peces,c2_peces) = NaN;
            end
%         end
    end % c1_peces
end % c2_peces

figure; 
ff = 10;
plot(delay,squeeze(correlmedia(ff,:,:))'); 
hold on; 
plot(delay_rel(ff,:)*framerate,maxcorr(ff,:),'o')

maxcorr=squeeze(maxcorr)+permute(squeeze(maxcorr),[2 1]);
delay_rel=(squeeze(delay_rel)-permute(squeeze(delay_rel),[2 1]));
delay_rel(maxcorr<delay_correlation_lim)=NaN;
delay_rel_abs=sign(delay_rel);
delay_rel_abs(isnan(delay_rel))=NaN;
means=nanmean(delay_rel,1)';
[~, orden]=sort(means,1,'descend'); % The last fish is #1, the first ('leader') #no_fish.
orden_out = NaN(no_fish,no_fish);
for ff=1:no_fish
    idx = find(orden==ff);
    orden_out(ff,ff) = idx;
end
correlmedia=mat2cell(correlmedia,ones(1,no_fish),ones(1,no_fish),ones(1,2*max_delay_frames+1));
%%
fstack=dbstack(3);
if ~isempty(fstack)
    funcinfo.Function = mfilename;
    funcinfo.callerFunction = fstack(1).name;
    funcinfo.no_fish = no_fish;
    funcinfo.XTick = edges(1:end-1);
    funcinfo.XTickLabel = cellfun(@(x) strtrim(x),cellstr((num2str(edges(1:end-1)')))','UniformOutput',false);
end
