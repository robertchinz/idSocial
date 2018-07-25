function [tr_out,vel,acc] = idSocial_auxiliaries_spline2arrays(tr,fnder_flag)
no_fish = numel(tr);

if nargin < 2 || isempty(fnder_flag)
    fnder_flag = false;
end

if isfield(tr{1},'indices_xx')
    no_frames = numel(tr{1}.indices_xx);
else
    no_frames = size(tr{1}.X,1);
end

tr_out = NaN(no_frames,no_fish,2);
vel = NaN(no_frames,no_fish,2);
acc = NaN(no_frames,no_fish,2);

% Calculate vel and acc directly from trajectories
% (values) instead of splines.


for ff=1:no_fish
    
    if isfield(tr{ff},'xx')
        trxtemp = fnval(tr{ff}.X,tr{ff}.xx);
        tr_out(tr{ff}.indices_xx,ff,1) = trxtemp;
    else
        tr_out(:,ff,1) = squeeze(tr{ff}.X);
    end
    if isfield(tr{ff},'yy')
        trytemp = fnval(tr{ff}.Y,tr{ff}.yy);
        tr_out(tr{ff}.indices_yy,ff,2) = trytemp;
    else
        tr_out(:,ff,2) = squeeze(tr{ff}.Y);
    end
    
end
if ~isfield(tr{ff},'velx') || ~isfield(tr{ff},'vely') || ...
        ~isfield(tr{ff},'accx') || ~isfield(tr{ff},'accy')
%     keyboard
    if ~fnder_flag
        vel(1:no_frames-1,:,:)=diff(tr_out,1,1);
        acc(1:no_frames-2,:,:)=diff(tr_out,2,1);
    else
        for ff=1:no_fish
            
            spvX = fnder(tr{ff}.X,1);
            vxtemp = fnval(spvX,tr{ff}.xx);
            vel(tr{ff}.indices_xx,ff,1) = vxtemp;
            spvY = fnder(tr{ff}.Y,1);
            vytemp = fnval(spvY,tr{ff}.yy);
            vel(tr{ff}.indices_yy,ff,2) = vytemp;
            
            spaX = fnder(tr{ff}.X,2);
            axtemp = fnval(spaX,tr{ff}.xx);
            acc(tr{ff}.indices_xx,ff,1) = axtemp;
            spaY = fnder(tr{ff}.Y,2);
            aytemp = fnval(spaY,tr{ff}.yy);
            acc(tr{ff}.indices_yy,ff,2) = aytemp;
        end
    end
else
    for ff=1:no_fish
        % Or: If vel and acc already exist:
        if isfield(tr{ff},'velx')
            vel(:,ff,1) = tr{ff}.velx;
        else
            %         spvX = fnder(tr{ff}.X,1);
            %         vxtemp = fnval(spvX,tr{ff}.xx);
            %         vel(tr{ff}.indices_xx,ff,1) = vxtemp;
        end
        
        if isfield(tr{ff},'vely')
            vel(:,ff,2) = tr{ff}.vely;
        else
            %         spvY = fnder(tr{ff}.Y,1);
            %         vytemp = fnval(spvY,tr{ff}.yy);
            %         vel(tr{ff}.indices_yy,ff,2) = vytemp;
        end
        
        if isfield(tr{ff},'accx')
            acc(:,ff,1) = tr{ff}.accx;
        else
            %         spaX = fnder(tr{ff}.X,2);
            %         axtemp = fnval(spaX,tr{ff}.xx);
            %         acc(tr{ff}.indices_xx,ff,1) = axtemp;
        end
        if isfield(tr{ff},'accy')
            acc(:,ff,2) = tr{ff}.accy;
        else
            %         spaY = fnder(tr{ff}.X,2);
            %         aytemp = fnval(spaY,tr{ff}.yy);
            %         acc(tr{ff}.indices_yy,ff,2) = aytemp;
        end
        
    end
end