function idSocial_auxiliaries_saveworkspace


thispath=mfilename('fullpath');
breakidx=strfind(thispath,'\');
savepath=[thispath(1:breakidx(end-1)) 'figures\'];

evalin('base',['save(''' savepath 'idSocial_workspace.mat'',''-v7.3'')']);
disp(['Saved workspace to ' savepath 'idSocial_workspace.mat'])