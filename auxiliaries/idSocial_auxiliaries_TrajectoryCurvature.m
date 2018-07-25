function curv = idSocial_auxiliaries_TrajectoryCurvature(trajectories)



if ~isempty(trajectories)
    [tr,~,~,no_frames,no_fish,no_dims] =idSocial_auxiliaries_formatInputTrajectory(trajectories,false);
    
    trajectories = NaN(no_frames,no_fish,no_dims);
    for ff=1:no_fish
        trajectories(:,ff,:) = tr(:,ff,ff,:);
    end
% (Slight) adaptation from Roger Stafford,
% http://www.mathworks.com/matlabcentral/answers/57194-how-to-find-the-sharp-turn-in-a-2d-line-curve#answer_69185

x1 = NaN(no_frames,no_fish,1);
x3 = NaN(no_frames,no_fish,1);
y1 = NaN(no_frames,no_fish,1);
y3 = NaN(no_frames,no_fish,1);

x1(2:end,:) = squeeze(trajectories(1:end-1,:,1));
x2 = squeeze(trajectories(1:end,:,1));
x3(1:end-1,:) = squeeze(trajectories(2:end,:,1));

y1(2:end,:) = squeeze(trajectories(1:end-1,:,2));
y2 = squeeze(trajectories(1:end,:,2));
y3(1:end-1,:) = squeeze(trajectories(2:end,:,2));

curv = 2*abs((x2-x1).*(y3-y1)-(x3-x1).*(y2-y1)) ./ ...
   sqrt(((x2-x1).^2+(y2-y1).^2).*((x3-x1).^2+(y3-y1).^2).*((x3-x2).^2+(y3-y2).^2));
else
    curv = NaN;
end