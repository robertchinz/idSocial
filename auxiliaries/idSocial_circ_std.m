function mean_std = idSocial_circ_std(ang)

vec=[cos(ang) sin(ang)];
R = nansum(vec,1);
R = R/sqrt(sum(R.^2,2));
mean_ang = atan2(R(2),R(1));
%%
% figure; plot(vec(:,1), vec(:,2),'.'); axis equal
% hold on
% plot([0 R(1)],[0,R(2)],'-r')