% 04-Jun-2012 16:16:14 Hago que se simetrice, para no tener sólo las
% correlaciones en la parte triangular superior.
% APE 28 may 12

function [correl,autocorr,delayed_dist]=trayectorias2correlacion(trayectorias,delay)
% Calculates correlation of movement directions; delay<0: focal follows neighbor; 
% delay>=0: focal leads. 

if nargin<2 || isempty(delay)
    delay=0;
end
no_dim = size(trayectorias,3);
vel=diff(trayectorias,1,1);
normavel=sqrt(sum(vel.^2,3));
velnorm=vel./repmat(normavel,[1 1 no_dim]);
n_peces=size(trayectorias,2);
n_frames=size(vel,1);
correl=NaN(n_frames,n_peces,n_peces);
delayed_dist=NaN(n_frames,n_peces,n_peces);
for c1_peces=1:n_peces
    for c2_peces=c1_peces+1:n_peces
        if delay>=0
            correl(:,c1_peces,c2_peces)=sum(velnorm(:,c1_peces,:).*[NaN(delay,1,no_dim) ; velnorm(1:end-delay,c2_peces,:)],3);
        else
            correl(:,c1_peces,c2_peces)=sum([NaN(-delay,1,no_dim) ; velnorm(1:end+delay,c1_peces,:)].* velnorm(1:end,c2_peces,:),3);
        end        
        correl(:,c2_peces,c1_peces)=correl(:,c1_peces,c2_peces);
        delayed_dist(:,c1_peces,c2_peces) = ...
            sqrt((trayectorias(1:end-1,c1_peces,1) - trayectorias(1:end-1,c2_peces,1)).^2 + (trayectorias(1:end-1,c1_peces,2) - trayectorias(1:end-1,c2_peces,2)).^2);
    end % c2_peces    
    correl(:,c1_peces,c1_peces)=1;
end % c1_peces
if nargout>=2
    autocorr=sum(velnorm(1:end-1,:,:).*velnorm(2:end,:,:),3);
end

