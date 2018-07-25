function outp=idSocial_velocityDensity(x,y,pos,vel,v1)


edg=0:max(vel(:))/(size(vel,1)/10):max(vel(:)); %v1=vel(fr-1,ff);
vbin=find(histc(v1,edg)); 
velfilter=vel(vel(:,1)>edg(vbin)& vel(:,1)<edg(vbin+1),1);
[f,xi] = ksdensity(velfilter);
outp=f([sqrt((x-pos(1))^2+(y-pos(2))^2)>xi(1:end-1) & ...
     sqrt((x-pos(1))^2+(y-pos(2))^2)<xi(2:end) false]);
if isempty(outp);
    outp=0;
end