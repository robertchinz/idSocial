function [def_options,def_types] = idSocial_auxiliaries_createDefOptions(description_flag)

if nargin<1 || isempty(description_flag)
    description_flag = false;
end

def_options.timeintervals_in_min= inf;
def_options.filter_focal_speedlimits_bl_per_s=[0 inf];
def_options.filter_neighbor_speedlimits_bl_per_s=[0 inf];
def_options.filter_neighbor_accelerationlimits_bl_per_s2 = [0 inf];
def_options.filter_focal_accelerationlimits_bl_per_s2  = [0 inf];
def_options.filter_distancelimits_bl=[0 inf];
def_options.filter_focal_rectangularROI = [-inf inf -inf inf];
def_options.filter_neighbor_rectangularROI = [-inf inf -inf inf];
% def_options.filter_focal_circularROI = [0 0 0 inf];
def_options.filter_focal_circularROI = {-Inf Inf {'BL';'Arena'}};
def_options.filter_neighbor_circularROI = {-Inf Inf {'BL';'Arena'}};
def_options.filter_AllMembersPresent = 0; 
def_options.filter_AllMembersPresentWorstFocal = 0; 
def_options.filter_AllMembersPresentWorstFocal = 0; 
def_options.filter_WorstIndividual = 0;
% def_options.filter_CustomFrames = [];
def_options.filter_framesWithAllNeighborsOnly = false;
def_options.filter_focal_spatial_sectors = {'all','lateral','frontal','behind'};
def_options.order_neighbors = false;
def_options.order_com = false;
def_options.random_data=true;
% def_options.temp_savepath = '<path>';
% def_options.parallel_processing = true;
% def_options.significance_between_groups = true;

if description_flag
    def_options(2).timeintervals_in_min='Intervals in minute used for binning data.';
    def_options(2).filter_focal_speedlimits_bl_per_s='Filter';
    def_options(2).filter_neighbor_speedlimits_bl_per_s='Filter';
    def_options(2).filter_neighbor_accelerationlimits_bl_per_s2 = 'Filter';
    def_options(2).filter_focal_accelerationlimits_bl_per_s2  = 'Filter';
    def_options(2).filter_distancelimits_bl='Filter';
    def_options(2).filter_focal_rectangularROI = 'Filter';
    def_options(2).filter_neighbor_rectangularROI = 'Filter';
    def_options(2).filter_focal_circularROI = 'Filter';
    def_options(2).filter_neighbor_circularROI = 'Filter';
    def_options(2).filter_AllMembersPresent = 'Filter';
    def_options(2).filter_AllMembersPresentWorstFocal = 'Filter';
    def_options(2).filter_AllMembersPresentWorstFocal = 'Filter';
    def_options(2).filter_WorstIndividual = 'Filter';
%     def_options(2).filter_CustomFrames = 'Filter';
    def_options(2).filter_framesWithAllNeighborsOnly = 'Filter';
    def_options(2).filter_focal_spatial_sectors = 'Filter';
    def_options(2).order_neighbors = 'Sort neighbors';
    def_options(2).order_com = 'Sort neighbors';
    def_options(2).random_data= 'Random data';
%     def_options(2).temp_savepath = 'Path';
%     def_options(2).parallel_processing = 'Activate parallel processing';
%     def_options(2).significance_between_groups = '';
end
%% Types
def_types.timeintervals_in_min={'numeric' ''};
def_types.filter_focal_speedlimits_bl_per_s={'numeric' ''};
def_types.filter_neighbor_speedlimits_bl_per_s={'numeric' ''};
def_types.filter_neighbor_accelerationlimits_bl_per_s2 ={ 'numeric' ''};
def_types.filter_focal_accelerationlimits_bl_per_s2  ={ 'numeric' ''};
def_types.filter_distancelimits_bl={'numeric' ''};
def_types.filter_focal_rectangularROI ={ 'numeric' ''};
def_types.filter_neighbor_rectangularROI ={ 'numeric' ''};
def_types.filter_focal_circularROI ={ 'numeric' ''};
def_types.filter_neighbor_circularROI ={ 'numeric' ''}; 
def_types.filter_AllMembersPresent ={ 'numeric' ''}; 
def_types.filter_AllMembersPresentWorstFocal ={ 'numeric' ''}; 
def_types.filter_AllMembersPresentWorstFocal ={ 'numeric' ''}; 
def_types.filter_WorstIndividual ={ 'numeric' ''};
% def_types.filter_CustomFrames ={ 'numeric' ''};
def_types.filter_framesWithAllNeighborsOnly = { 'numeric' ''};
def_types.filter_focal_spatial_sectors ={ 'numeric' ''};
def_types.order_neighbors ={ 'logical' ''};
def_types.order_com = {'logical' ''};
def_types.random_data={'logical' ''};
% def_types.temp_savepath ={ '' ''}; 
% def_types.parallel_processing ={ 'logical' ''};
% def_types.significance_between_groups ={ 'logical' ''};