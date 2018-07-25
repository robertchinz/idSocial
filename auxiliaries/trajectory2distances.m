function distances=trajectory2distances(trayectorias)

no_frames=size(trayectorias,1);
no_fish=size(trayectorias,2);
distances=NaN(no_fish,no_fish,no_frames);
for ffish=1:no_fish
    for nfish=[1:ffish-1 ffish+1:no_fish]
        distances(ffish,nfish,:)=sqrt(sum((trayectorias(:,ffish,:)-trayectorias(:,nfish,:)).^2,3));
    end 
end