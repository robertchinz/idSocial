% trayectorias must be a n_frames x n_individuals x 2 matrix.

% 21-Oct-2011 12:09:52 Corrijo, por si hay NaN al principio
% APE 20 oct 11

function trayectorias=trayectorias2trayectorias_interp(trayectorias)

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