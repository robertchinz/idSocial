function [hi,prob] = idSocial_auxiliaries_kerneldensity(val,edges_kds,kds_bandwidth,method)

if nargin<4 || isempty(method)
    method='ksdensity_epanechnikov';
end
if size(val,2)>1 && size(val,1)==1
    val=val';
end

cx=repmat(val,[1 numel(edges_kds)]);
edgesX=repmat(edges_kds,[size(val,1) 1]);
% keyboard
if ~strcmpi(method,'hist')
    if strcmpi(method,'ksdensity_epanechnikov')
        if nargout == 1 
            % Note: The bandwidth parameter bw in ksdensity
            % gives the same results as 2*bw when calculating it myself IN CASE OF EPANECHNIKOV.
            % For a normal/gaussian kernel it is the same.
            hi = ksdensity(val,edges_kds,'kernel','epanechnikov','width',.5*kds_bandwidth);
            prob = 1;
        else
            prob=.75*(1-(1/(kds_bandwidth)*(cx-edgesX)).^2)/(kds_bandwidth); % Epanechnikov kernel
            prob(prob<0)=NaN;
            norm = nansum(prob,2);
            norm(norm==0) = NaN;
            prob = prob./repmat(norm,[1 size(prob,2)]);
            hi=nansum(prob,1);
%             figure; plot(edges_kds,hi1,'+'); hold on; plot(edges_kds,hi2/sum(hi2)/(edges_kds(2)-edges_kds(1)),'r*')
        end
    elseif strcmpi(method,'ksdensity_triangular')
        if nargout == 1 
            hi = ksdensity(val,edges_kds,'kernel','triangle','width',kds_bandwidth);
            prob = 1;
        else
            prob=(1-1/kds_bandwidth*abs(cx-edgesX))/kds_bandwidth; % Triangular kernel
            prob(prob<0)=NaN;
            hi=nansum(prob,1);
        end
    elseif strcmpi(method,'ksdensity_gauss')
        if nargout == 1 
            hi = ksdensity(val,edges_kds,'kernel','normal','width',kds_bandwidth);
            prob = 1;
        else
            prob=exp(...
                -.5*((cx-edgesX)/kds_bandwidth).^2 )/(sqrt(2*pi)*kds_bandwidth);
            hi=nansum(prob,1);
        end
    end
    
%     hi=hi/sum(hi);
else
    hi = histc(val,edges_kds);
    hi = hi/sum(hi);
    prob = [];
end