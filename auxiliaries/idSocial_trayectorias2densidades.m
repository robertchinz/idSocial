% sólo metemos un frame para que saque las densidades

function mapa_densidades=idSocial_trayectorias2densidades(trayectorias_pixels,R_influencia,Lx,Ly,idpez,sello)

% R_influencia=500;
if nargin<6 || isempty(sello)
    sello=pos2sello(R_influencia);
end

Nindiv=size(trayectorias_pixels,2);
% Nit=size(trayectorias,1);

% xmin=min(min(trayectorias(:,:,1)));
% xmax=max(max(trayectorias(:,:,1)));
% ymin=min(min(trayectorias(:,:,2)));
% ymax=max(max(trayectorias(:,:,2)));

% trayectorias_origen(:,:,1)=trayectorias(:,:,1)-xmin;
% trayectorias_origen(:,:,2)=trayectorias(:,:,2)-ymin;
% 
% trayectorias_pixels=round(trayectorias_origen);
% 
% Lx=ceil(xmax)-floor(xmin);
% Ly=ceil(ymax)-floor(ymin);

% figure
t=1;
% for t=1:Nit
    mapa_densidades=zeros(Ly,Lx,'uint8');
%     clf
   
    
    
    for i=1:Nindiv
        
        if i~= idpez   %idpez
            
%             disp('caca que has cambiado para pintar densidades!!!')
        
        inicio_x=trayectorias_pixels(t,i,1)-R_influencia;
        fin_x=trayectorias_pixels(t,i,1)+R_influencia;
        inicio_y=trayectorias_pixels(t,i,2)-R_influencia;
        fin_y=trayectorias_pixels(t,i,2)+R_influencia;

        
        
        sello_it=sello;
        if trayectorias_pixels(t,i,1)+R_influencia>Lx
            sello_it(:,R_influencia+1+Lx-trayectorias_pixels(t,i,1)+1:2*R_influencia+1)=[];
            fin_x=Lx;
        end
        
        
        if trayectorias_pixels(t,i,1)-R_influencia<1        
            sello_it(:,1:R_influencia+1-trayectorias_pixels(t,i,1))=[];
            inicio_x=1;
        end

        if trayectorias_pixels(t,i,2)+R_influencia>Ly
            sello_it(R_influencia+1+Ly-trayectorias_pixels(t,i,2)+1:2*R_influencia+1,:)=[];
            fin_y=Ly;
        end
        
        if trayectorias_pixels(t,i,2)-R_influencia<1        
            sello_it(1:R_influencia+1-trayectorias_pixels(t,i,2),:)=[];
            inicio_y=1;
        end
        
        capa_indiv=zeros(Ly,Lx,'uint8');
%         try
        capa_indiv(inicio_y:fin_y,inicio_x:fin_x)=sello_it;
%         catch
%             keyboard
%         end
        
        mapa_densidades=mapa_densidades+capa_indiv;
         
        end
        
        
        
        
        
    end   

% mapa_social=mapa_densidades~=0;
% max_densidad=max(max(mapa_densidades));
    


% posicion_Npeces=zeros(Ly,Lx,Nindiv-1);
% 
% for i=1:Nindiv-1
%     posicion_Npeces(:,:,i)=(mapa_densidades==i);
% end



%     imagesc(mapa_densidades)
%     axis equal
%     drawnow
%     keyboard
% end





