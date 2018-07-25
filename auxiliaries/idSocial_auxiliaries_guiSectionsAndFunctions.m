function [sectionNames funcCell]=idSocial_auxiliaries_guiSectionsAndFunctions(dirname) 

if nargin<1 || isempty(dirname)
    dr=textscan(mfilename('fullpath'),'%s','Delimiter',[filesep]);
    dr2=cellfun(@(x) [x filesep],dr{1},'UniformOutput',false);
    dirname=[dr2{1:end-2}];
end
    
d = dir(dirname);
isub = [d(:).isdir]; %# returns logical vector
nameFolds = {d(isub).name}';
nameFolds=nameFolds(cellfun(@(x) strcmp(x(1),'~'),nameFolds));
sectionNames=cellfun(@(x) regexprep(x(2:end),'(\<[a-z])','${upper($1)}'),nameFolds,'UniformOutput',false);

funcCell=cell(size(nameFolds,1),3);
for k=1:size(nameFolds,1)
    subdir=[dirname nameFolds{k} filesep];
    sd = dir(subdir);
    isub = [sd(:).isdir];
    funcs = {sd(~isub).name};
    functs=cellfun(@(x) textscan(x,'%s','Delimiter','._'),funcs);
    ismfile=cellfun(@(x) strcmp(x{end},'m'),functs);
    isidSocial=cellfun(@(x) strcmpi(x{1},'idSocial'),functs);
    isInSection=cellfun(@(x) strcmpi(x{2},sectionNames{k}),functs);
    if ~isempty(funcs(ismfile&isidSocial&isInSection)')
        funcCell{k,1}=funcs(ismfile&isidSocial&isInSection)';
        funcNames=cellfun(@(x) textscan(x,'%s','Delimiter','_.'),funcCell{k,1});
        funcNames=cellfun(@(x)regexprep([x{3:end-1}],'(\<[a-z])','${upper($1)}'),funcNames,'UniformOutput',false);
        funcCell{k,2}=funcNames;
    end
    for k2=1:size(funcCell{k,1},1)
        [~, funcCell{k,3}{k2}]=idSocial_auxiliaries_guiDefOptions( [subdir funcCell{k,1}{k2}]);
    end
end
emptySec=cellfun(@(x) isempty(x), funcCell(:,1));
funcCell(emptySec,:)=[];
sectionNames(emptySec,:)=[];
