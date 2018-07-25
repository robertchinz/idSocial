function varargout=idSocialOptionsEditor(options)

% out = [];
varargout = cell(1,1);
% Add Java library to dynamic Java classpath

% javaaddpath(['C:\Users\Robert\workspace\ExampleWindow' '\ExampleWindow.jar']);
% javaaddpath(['C:\Users\Robert\workspace\ExampleWindow' '\ExampleWindow.jar']);
javaaddpath(['C:\Users\Robert\workspace\idSocialOptionsEditorUI.jar']);
javaaddpath(['C:\Users\Robert\workspace\miglayout-4.0-swing.jar']);

% gui=guidata(fh);

% javaaddpath C:\Users\Robert\workspace\ExampleWindowMatlab1_7\bin\examplewindow

% Get example Java window from the library
jFrame = idSocialOptionsEditor.IdSocialOptsEdit();

% Display the Java window
jFrame.setVisible(true);

% keyboard
% OK Button
showOkButton = handle(jFrame.getOkButton(), 'CallbackProperties');
showCancelButton = handle(jFrame.getCancelButton(), 'CallbackProperties');


set(showOkButton, 'ActionPerformedCallback', @(src,evnt) okCallback(src,evnt));
set(showCancelButton, 'ActionPerformedCallback', @(src,evnt) cancelCallback(src,evnt));

% Get option strings and values
optString = fieldnames(options);
optVals = struct2cell(options);

jFrame.setTableData(optString,optVals);

    function okCallback(src,evnt)
        newData = jFrame.getTableData(1);
        varargout{1}  = newData;
        
        closeWindow;
    end

    function cancelCallback(src,evnt)
        jFrame.dispose
    end
    function closeWindow
        
        jFrame.dispose
    end

if isdeployed
    waitfor(jFrame);
end
end
% function okCallback(src,evnt,fh,jFrame)
% gui=guidata(fh);
% newData = jFrame.getTableData(1);
% gui.newOptions = newData;
% guidata(fh,gui);
% jFrame.dispose
% end
%
% function cancelCallback(src,evnt,jFrame)
% jFrame.dispose
% end
%%
% jFrame.dispose();