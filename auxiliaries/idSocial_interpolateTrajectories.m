function trayectorias=idSocial_interpolateTrajectories(trayectorias,mode,maxlen)

if nargin<2 || isempty(mode)
    mode='linear';
% else
%     mode='spline';
end
if nargin<3 || isempty(maxlen)
    maxlen=inf;
end

no_fish=size(trayectorias,2);

nanidx=cell(no_fish,1);
nanstartidx=cell(no_fish,1);
nanendidx=cell(no_fish,1);
nanlength=cell(no_fish,1);
for ff=1:no_fish
    nanidx{ff}=find(isnan(trayectorias(:,ff,1)));
    if ~all(any(isnan(trayectorias(:,ff,:)),3),1) && ~isempty(nanidx{ff}) && length(nanidx{ff})>1 
        nanstartidx{ff}=nanidx{ff}([1 find(diff(nanidx{ff})>1)'+1]);
        nanendidx{ff}=nanidx{ff}(find(diff(nanidx{ff})>1)');
        if ~isempty(nanendidx{ff})
            try
                nanlength{ff}=nanendidx{ff}-nanstartidx{ff}(1:end-1)+1;
            catch
                keyboard
            end
            
            nanstartidx{ff}= nanstartidx{ff}(nanlength{ff}>maxlen);
            nanendidx{ff}=nanendidx{ff}(nanlength{ff}>maxlen);
            %         nanidx{ff}=isnan(trayectorias(:,ff,1));
            %         nanstartidx{ff}=[0 diff(nanidx{ff})'>0];
            %         nanendidx{ff}=diff(nanidx{ff})'<0;
            %     nanendidx{ff}=nanidx{ff}(find(diff(nanidx{ff})>1)');
            %     nanlength{ff}=find(nanendidx{ff})-find(nanstartidx{ff}(1:end-1))+1;
            for k=1:length(nanstartidx{ff})
                trayectorias(nanstartidx{ff}(k):nanendidx{ff}(k),ff,:)=inf;
            end
        else
            nanendidx{ff}=nanstartidx{ff};
        end
    end
end


if strcmpi(mode,'spline')
    
    trayectorias(isnan(trayectorias)) = ...
        interp1(find(~isnan(trayectorias)), trayectorias(~isnan(trayectorias)), find(isnan(trayectorias)),'spline');
    
else
    % The following code is taken from the function
    % trayectorias2trayectorias_interp(trayectorias) by Alfonso Pérez Escudero
    % trayectorias must be a n_frames x n_individuals x 2 matrix.
    
    % 21-Oct-2011 12:09:52 Corrijo, por si hay NaN al principio
    % APE 20 oct 11
    
    %
    
    n_frames=size(trayectorias,1);
    n_peces=size(trayectorias,2);
    for c_peces=1:n_peces
        c_frames=1;
        while c_frames<n_frames
            % Avanza hacia el primer hueco
            while c_frames<=n_frames && ~isnan(trayectorias(c_frames,c_peces))
                c_frames=c_frames+1;
            end
            if c_frames<=n_frames
                frame1=c_frames;
                if c_frames>1
                    pos0=squeeze(trayectorias(c_frames-1,c_peces,:))';
                else
                    pos0=[];
                end
                %     v0=squeeze(diff(trayectorias(c_frames-2:c_frames-1,c_peces,:),1,1))';
                n_huecos=0;
                % Avanza hasta el final del hueco
                while c_frames<=n_frames && isnan(trayectorias(c_frames,c_peces))
                    n_huecos=n_huecos+1;
                    c_frames=c_frames+1;
                end
                if ~isempty(pos0) && c_frames<=n_frames
                    frame2=c_frames-1;
                    posf=squeeze(trayectorias(c_frames,c_peces,:))';
                    %     vf=squeeze(diff(trayectorias(c_frames:c_frames+1,c_peces,:),1,1))';
                    % Interpolación
                    pasos=0:1/(n_huecos+1):1;
                    pasos=pasos(2:end-1);
                    %     pos=NaN(n_huecos,2);
                    trayectorias(frame1:frame2,c_peces,1)=pos0(1)+pasos*(posf(1)-pos0(1));
                    trayectorias(frame1:frame2,c_peces,2)=pos0(2)+pasos*(posf(2)-pos0(2));
                end
            end
        end % while quedan frames
    end % c_peces
end
for ff=1:no_fish
    for k=1:length(nanstartidx{ff})
        trayectorias(nanstartidx{ff}(k):nanendidx{ff}(k),ff,:)=NaN;
    end
end