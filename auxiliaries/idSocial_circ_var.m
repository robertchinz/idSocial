function circvar = idSocial_circ_var(ang)
% See A MATLAB Toolbox for Circular Statistics, Philipp Berens

vec=[cos(ang) sin(ang)];
R = nanmean(vec,1);
circvar = 1-sqrt(sum(R.^2,2));

%%
% figure; plot(vec(:,1), vec(:,2),'.'); axis equal
% hold on
% plot([0 R(1)],[0,R(2)],'-r')