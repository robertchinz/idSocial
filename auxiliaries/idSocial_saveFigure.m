function idSocial_saveFigure(filename,project_path)

if nargin<2 || isempty(project_path)
    project_path=pwd;
end

thispath=project_path;
if ~strcmp(thispath(end),'\'); thispath=[thispath '\']; end;

visible=get(gcf,'visible');
set(gcf,'Color','w')

screen_size = get(0, 'ScreenSize');

% pause(.3)
% jFrame = get(handle(gcf),'JavaFrame');
% jFrame.setMaximized(true)
if strcmpi(visible,'Off');
    %     pos=get(gcf,'Position');
    set(gcf,'Position',[-screen_size(3)-10 -screen_size(4)-10 screen_size(3) screen_size(4)]);
    %     pause(.3)
    %     jFrame = get(handle(gcf),'JavaFrame');
    %     jFrame.setMaximized(true)
    set(gcf,'Visible','On')
    snapnow
    set(gcf,'Visible','Off')
end

if ~isempty(project_path)
    titlestring=filename;
    %     breakidx=strfind(thispath,'\');
    %     savepath=[thispath(1:breakidx(end-1)) 'figures\'];
    
    savepath=[thispath 'figures\'];
    
    if exist(savepath,'dir')==0
        try
        mkdir(savepath)
        disp(['Created ' savepath])
        catch
            disp(['Could not create ' savepath]) 
            
        end
    end
    %     savepath_latex=[thispath(1:breakidx(end-1)) 'figures_latex\'];
    savepath_latex=[thispath 'figures_latex\'];
    
    if exist(savepath_latex,'dir')==0
        try
        mkdir(savepath_latex)
        disp(['Created ' savepath_latex])
        catch
            disp(['Could not create ' savepath]) 
        end
    end
    
    
    datestring=datestr(now, 30);
    
    if strcmpi(visible,'Off');
        set(gcf,'Position',[-screen_size(3)-10 -screen_size(4)-10 screen_size(3) screen_size(4)]);
        
        set(gcf,'Visible','On')
        
    elseif strcmpi(visible,'On');
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
        pause(.3)
        jFrame = get(handle(gcf),'JavaFrame');
        jFrame.setMaximized(true)
    end
    
    pause(.3)
    set(gcf,'PaperPositionMode','auto')
    try
        if nargin<1 || isempty(filename)
            export_fig([savepath_latex '_' datestring],'-transparent','-pdf','-png')
        else
            export_fig([savepath_latex titlestring],'-transparent','-pdf','-png')
        end
    catch
        disp(['FAILED TO SAVE ' savepath_latex titlestring])
    end
    try
        export_fig([savepath titlestring '_' datestring],'-transparent','-png')
        saveas(gcf, [savepath titlestring '_' datestring],'fig');
        disp(['Saved ' savepath titlestring '_' datestring '.pdf, ~.png and ~.fig' ]);
    catch
        disp(['FAILED TO SAVE ' savepath titlestring '_' datestring])
    end
    set(gcf,'Visible',visible)
    %     try
    %     pause(.6)
    %     saveas(gcf, [savepath titlestring '_' datestring],'fig');
    %      catch
    %         keyboard
    %     end
    
    % saveas(gcf, [savepath titlestring],'fig');
    
    
    % set(gcf,'PaperUnits','normalized');
    % % set(gcf,'PaperPosition', [0 0 1 1]);
    % set(gcf, 'PaperPositionMode','auto')
    % print(gcf,[savepath titlestring '_' datestring],'-djpeg','-r600');
    % print(gcf,[savepath titlestring],'-dpdf','-r600');
    
    
end
% if strcmpi(visible,'Off')
%     close(gcf) % It only saves the images, does not display them.
% else
%     set(gcf,'Visible','On')
% end

