function [trajectory_location, directory, trloc_temp]= idSocial_recursiveFolderSearch(directory,fname,trajectory_location)

if nargin < 3 || isempty(trajectory_location)
    trajectory_location = cell(1,1);
end

files = dir(directory);
fnames = vertcat({files.name});
dirs = fnames(vertcat(files.isdir)' & ~strcmp(fnames,'.') & ~strcmp(fnames,'..'));

fnameidx = strcmpi(fnames,fname);

if strcmp(directory(end),'\')
    directory = directory(1:end-1);
end

if any(fnameidx) && isempty(dirs)
   trajectory_location = [directory '\' fname];
   return;
elseif ~isempty(dirs)
    trajectory_location = cell(numel(dirs),1);
    for k=1:numel(dirs)
        directory_new = [directory '\' dirs{k}];
        trajectory_location{k} = idSocial_recursiveFolderSearch(directory_new,fname,trajectory_location{k});
        
    end
end



