function tray_rand=idSocial_auxiliaries_randomizeTrajectory(trayectorias,no_randomizations)

no_fish=size(trayectorias,2);
if nargin<2 || isempty(no_randomizations)
    no_randomizations=1;
end
tray_rand=NaN(size(trayectorias,1),size(trayectorias,2)*no_randomizations,2);
for rep=1:no_randomizations
    for c=1:no_fish
        tray_rand(:,(rep-1)*no_fish+c,:)=tray_rand(randperm(size(trayectorias,1)),c,:);
    end
end
