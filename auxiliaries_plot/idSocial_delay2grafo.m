% 02-Aug-2012 10:08:22 Añado la posibilidad de que salga bonito (y sin
% ejes)
% 02-Aug-2012 10:01:50 Cambio la forma de calcular las posiciones.
% 01-Aug-2012 16:51:03 Mejoras de estética, y añado el umbral
% APE 29 may 12

% Seguro que hay una forma analítica de hacerlo invirtiendo una matriz.
% Pero lo voy a hacer cutre, calculando los centros de masas
% iterativamente.

function h=idSocial_delay2grafo(delay,superanumbral,bonito,hsignificant,colorines)
% DELAY is an N-by-N matrix with elements (i,j).
% Direction of arrows: (i,j)<0 (and (j,i)>0): Arrow from i 
% to j.  
% Vertical position of node i: Determined by 
% Mean([-(i,1),-(i,2),...,-(i,N)])   
% set(gcf,'Renderer','painters')
if nargin<2 || isempty(superanumbral)
    superanumbral=true(size(delay));
end
if nargin<4 || isempty(hsignificant)
    hsignificant=false(size(delay));
end
if nargin<3 || isempty(bonito)
    bonito=false;
end
if nargin<5 || isempty(colorines)
        colorines=[0 0 0; 0 1 0; 0 0 1; 1 1 0; 0 1 1; 1 0 1; .5 0 1; 1 0 0; 0 1 0; 0 0 1; 1 0 1; 1 .5 .5; 0 1 1; 1 1 0;  .25 .5 .5; 1 .5 0; 0 .6 0; .5 .2 1; .5 1 1; 1 1 0];
end
linewidth=1;
tipangle=10;
arrowlength=10;
radio=.09;
fontsize = 8;

n_bichos=size(delay,1);
% pos=[0 delay(1,2:end)]; % Pongo el primero como referencia
% pos_old=NaN(1,n_bichos);
% while ~all(abs(pos-pos_old)<10^-10)
% %     pos-pos_old
%     pos_old=pos;
%     pos_mat=repmat(pos,[n_bichos 1]);
%     pos_mat(1:(n_bichos+1):n_bichos^2)=0; % Anulo la diagonal
%     pos=(sum(pos_mat-delay,2)/(n_bichos-1))';
% end

delay_act=delay;
delay_act(~superanumbral)=NaN;
pos=nanmean(delay_act,1);
pos(isnan(pos))=nanmin(pos);
pos(isnan(pos))=0;

% x=rand(1,n_bichos);
x=0:1/(n_bichos-1):1;
% figure
if ~bonito % Versión cutre con ejes
    axis([-.05 1.05 min(pos)-.05 max(pos)+.05]);
    hold on
    for c1_bichos=1:n_bichos
        for c2_bichos=c1_bichos+1:n_bichos
            if superanumbral(c1_bichos,c2_bichos)
                if delay(c1_bichos,c2_bichos)<0
                    inicio=[x(c1_bichos) pos(c1_bichos)];
                    final=[x(c2_bichos) pos(c2_bichos)];
                else
                    inicio=[x(c2_bichos) pos(c2_bichos)];
                    final=[x(c1_bichos) pos(c1_bichos)];
                end
                arrow(inicio,final)
            end % if supera umbral
        end % c2_bichos
    end % c1_bichos
    
    for c1_bichos=1:n_bichos
        text(x(c1_bichos),pos(c1_bichos),num2str(c1_bichos),'FontSize',15,'FontWeight','bold','Color','r','HorizontalAlignment','center')
    end % c1_bichos
    
else % Versión bonita sin ejes
    %     colorines=[0 0 0 ; .5 0 1; 1 .5 0; .2 .7 .7; 0 .6 0; 0 0 1; .5 .5 0; 1 0 0; .25 .5 .5; .5 1 1; 1 1 0];
% [0 0 0 ; 1 0 0 ; 0 0 1;  .5 .5 0;  .25 .5 .5; 1 .5 0; 1 0 1; 0 .6 0; .5 .2 1; .5 1 1; 1 1 0];
    
    pos=pos-min(pos);
    if max(pos)>0; pos=pos/max(pos);end;
    axis([-radio 1+radio -radio 1+radio]);
    hold on    
    for c1_bichos=1:n_bichos
        for c2_bichos=c1_bichos+1:n_bichos
            if superanumbral(c1_bichos,c2_bichos)
                if delay(c1_bichos,c2_bichos)<0 && delay(c2_bichos,c1_bichos)>0
                    inicio=[x(c1_bichos) pos(c1_bichos)];
                    final=[x(c2_bichos) pos(c2_bichos)];
                    vector=final-inicio;
                    vector=vector/norm(vector);
                    if hsignificant(c1_bichos,c2_bichos)
                        colorflechas=[.2 .7 .2];
                    else
                        colorflechas=[.7 .7 .8];
                    end
                    arrow([inicio+radio*vector -1],[final-radio*vector -1],'LineWidth',linewidth,'EdgeColor',colorflechas,'FaceColor',colorflechas,'TipAngle',tipangle,'Length',arrowlength)
                    
                elseif delay(c1_bichos,c2_bichos)>0 && delay(c2_bichos,c1_bichos)<0
                    inicio=[x(c2_bichos) pos(c2_bichos)];
                    final=[x(c1_bichos) pos(c1_bichos)];
                    vector=final-inicio;
                    vector=vector/norm(vector);
                    if hsignificant(c1_bichos,c2_bichos)
                        colorflechas=[.2 .7 .2];
                    else
                        colorflechas=[.7 .7 .8];
                    end
                    arrow([inicio+radio*vector -1],[final-radio*vector -1],'LineWidth',linewidth,'EdgeColor',colorflechas,'FaceColor',colorflechas,'TipAngle',tipangle,'Length',arrowlength)
%                     uistack(ah,'bottom');
elseif delay(c1_bichos,c2_bichos)==delay(c2_bichos,c1_bichos)
                    inicio=[x(c2_bichos) pos(c2_bichos)];
                    final=[x(c1_bichos) pos(c1_bichos)];
                    vector=final-inicio;
                    vector=vector/norm(vector);
                    if hsignificant(c1_bichos,c2_bichos)
                        colorflechas=[.2 .7 .2];
                    else
                        colorflechas=[.7 .7 .8];
                    end
                    line([inicio(1)+radio*vector(1) final(1)-radio*vector(1) ] ,...
                        [inicio(2)+radio*vector(2) final(2)-radio*vector(2) ],'LineWidth',linewidth,'Color',colorflechas)
                    
                end
            end % if supera umbral
        end % c2_bichos
    end % c1_bichos
    
    for c1_bichos=1:n_bichos
          ra=rectangle('Position',[(x(c1_bichos)-radio) (pos(c1_bichos)-radio) 2*radio 2*radio],'Curvature',[1 1],'LineWidth',linewidth,'FaceColor',colorines(c1_bichos,:));
%            ra=plot((x(c1_bichos)), (pos(c1_bichos)),'.','MarkerSize',45,'LineWidth',linewidth,'Color',colorines(c1_bichos,:));
%                   ra=patch([x(c1_bichos)+radio*cos(0:.001:2*pi+.001)], ...
%                       [pos(c1_bichos)+radio*sin(0:.001:2*pi+.001)], ...
%                       [ones(1,numel(0:.001:2*pi+.001))*1],'FaceAlpha',1,'LineWidth',linewidth,'FaceColor',colorines(c1_bichos,:));
%                 uistack(ra,'top')
%         rectangle('Position',[max(x(c1_bichos)-radio,0) max(pos(c1_bichos)-radio,0) 2*radio 2*radio],'Curvature',[1 1],'LineWidth',linewidth,'FaceColor',colorines(c1_bichos+1,:))
        text((x(c1_bichos)),(pos(c1_bichos)),num2str(c1_bichos),'FontSize',fontsize,'FontWeight','bold','Color','k','HorizontalAlignment','center','VerticalAlignment','middle')
    end % c1_bichos
    set(gca,'DataAspectRatio',[1 1 1],'XTick',[],'YTick',[],'XColor','w','YColor','w')
end
h=gca;