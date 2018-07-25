function outCell = idSocial_recursiveDir(dr,filepattern,outCell)
if nargin<2 || isempty(filepattern)
    filepattern = '';
end
if nargin<3 || isempty(outCell)
    outCell = {};
end

flist = dir(dr);
fnames = {flist.name};
isDir = [flist.isdir];

onlyDots = cellfun(@(x) strcmp(x,'.')|strcmp(x,'..'),fnames);
fnames = fnames(~onlyDots);
isDir = isDir(~onlyDots);

if any(isDir)
    outCellTemp = cell(numel(fnames),1);
    for k = 1:numel(fnames)
        outCellTemp{k}=idSocial_recursiveDir([dr filesep fnames{k}],filepattern,outCell);
    end
    outCell = vertcat(outCell,vertcat(outCellTemp{:}));
else
    flist = dir([dr filesep filepattern]);
    fnames = {flist.name};
    onlyDots = cellfun(@(x) strcmp(x,'.')|strcmp(x,'..'),fnames);
    fnames = fnames(~onlyDots);
    for k = 1:numel(fnames)
        outCell=vertcat(outCell,([dr filesep fnames{k}]));
    end
    return;
end