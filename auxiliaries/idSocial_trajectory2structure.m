function [trajectory_list trajectory_cell]=idSocial_trajectory2structure(directory)

if nargin<1 || isempty(directory)
    
    [trajectory_list, directory]=uigetfile('*.mat','Select trajectory file','Multiselect','On');
    if ~iscell(trajectory_list)
        templist=trajectory_list;
        trajectory_list=cell(1,1);
         trajectory_list{1}= templist;
    end
    trajectory_list=trajectory_list';
    
else
    
    listing=dir([directory '\*tra*ector*.mat']);
    trajectory_list={listing.name}';
end

trajectory_cell=cell(size(trajectory_list,1),1);
for f=1:size(trajectory_list,1)
    load([directory '\' trajectory_list{f}]);
    if exist('trayectorias','var')==1
        trajectory_cell{f}=trayectorias;
    end
end


