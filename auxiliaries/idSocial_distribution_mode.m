function [mode_loc, hitemp]=idSocial_distribution_mode(vec,edges,bandwidth,method,location)

if nargin < 3 || isempty(bandwidth)
    bandwidth = 1;
end
if nargin < 4 || isempty(method)
    method = 'ksdensity_epanechnikov';
end
if nargin < 5 || isempty(location)
    location = [];
end

if ~isempty(vec) && ~all(isnan(vec))
    if size(vec,1) == 1 && size(vec,2)>1
            vec=vec';
    end
        [hitemp, prob] = idSocial_auxiliaries_kerneldensity(vec,edges,bandwidth,method);
        clear prob;
        if isempty(location)
%             [~,locs] = findpeaks(hitemp,'sortstr','descend','npeaks',1);
            [~,locs] = idSocial_auxiliaries_findpeaks(hitemp,'descend',1);
            mode_loc = edges(locs);
        elseif strcmpi(location,'first')
%             [~,locs] = findpeaks(hitemp);
            [~,locs] = idSocial_auxiliaries_findpeaks(hitemp);
            mode_loc = edges(min(locs));
        elseif strcmpi(location,'last')
%             [~,locs] = findpeaks(hitemp);
            [~,locs] = idSocial_auxiliaries_findpeaks(hitemp);
            mode_loc = edges(max(locs));
        end
        
else
    mode_loc = NaN;
end
% if ~isempty(vec) && ~all(isnan(vec))
%     if strcmpi(method,'ksdensity')
%         hitemp = idSocial_auxiliaries_kerneldensity(vec,edges_kds,bandwidth,method);
% %         hitemp = ksdensity(vec,edges,'width',bandwidth);
%         [~,locs] = findpeaks(hitemp,'sortstr','descend','npeaks',1);
%         mode_loc = edges(locs);
%     elseif strcmpi(method,'ksdensity_epanechnikov')
%         if size(vec,1) == 1 && size(vec,2)>1
%             vec=vec';
%         end
% %         edges_kds=edges(1:end-1);
% %         cx=repmat(vec,[1 numel(edges_kds)]);
% %         edgesX=repmat(edges_kds,[size(vec,1) 1]);
% % %         try
% %         prob=.75*(1-(cx-edgesX).^2); % Epanechnikov kernel
% % %         catch
% % %             keyboard
% % %         end
% %         prob(prob<0)=NaN;
% %         hitemp=nansum(prob,1);
% %         
%         hitemp = idSocial_auxiliaries_kerneldensity(vec,edges_kds,bandwidth,method);
%         
%         [~,locs] = findpeaks(hitemp,'sortstr','descend','npeaks',1);
%         mode_loc = edges(locs);
%     end
% else
%     mode_loc = NaN;
% end
