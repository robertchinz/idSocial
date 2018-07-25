function plot_mode
% Each idSocial function creates a cell array ORIG_CELL,
% which is saved in a file at the location given in
% data.('name of the actual function').output.
% 
% ORIG_CELL contains the results of all the calculations for
% each group, subset (e.g., day), trial (single trajectory),
% time step, focal, neighbor, and additional fields 
% depending on the actual function. Thus, OUTPUT is a cell
% array with dimensions
%
% (#Groups -by- #Subsets -by- #Trials -by- #Individuals 
%   -by- #Individuals -by- ...).
%
% NOTE: In each group there may be a different number of 
% subsets, and for each subset a different number of trials
% etc. Therefore, dimension sizes here refer to the maximum
% number of subsets in all groups, and analogously for 
% trials, individuals etc.
%
% In general, we want to plot and compare different pieces 
% of all the data. This is done by the setting PLOT_MODE.
% PLOT_MODE uses string specifiers to distinguish subsets of 
% the data:
% 'Group', 'Day', 'Trial', 'Time', 'Focal', 'Neighbor', 
%                                           'EdgeX', 'EdgeY'
%
% PLOT_MODE is a structure with the following fields: 
%       XAXIS:  Specifies the x-axis: 
%                   plot_mode.xaxis = {'Time'}
%               will plot the data versus time.
%               Some functions do not permit changing the 
%               axes, for example in case of distributions/
%               histograms, the x-axis is set to 'EdgeX' 
%               automatically, corresponding to the bin 
%               edges set before running the function (and,
%               analogously, the y-axis is set to 'EdgeY'
%               in case of bivariate histograms/maps).
%
%       YAXIS:  see XAXIS. Only applicable to bivariate 
%               histograms/maps.
%
%       FILTER: Specifies/filters data for calculations. 
%               Filter is a cell array with as many 
%               rows as the number of data subsets you want 
%               to compare (e.g., the number of legend 
%               entries in a plot). Each row has the form
%
%                    {'Specifier 1' [indices] ...
%                          'Specifier 2' [indices] ...}
%               where 'Specifier' is a string specifier 
%               from the list above, and 'indices' determine
%               the data indices which will be used for the
%               calculations. E.g., 
%
%                    plot_mode.filter =
%                           {'Group' 1 'Trial' 1; ...
%                           'Group' 2 'Trial' 3 }
%
%               corresponds to data from Group 1, 
%               Trial 1 in the first plot, and Group 2, 
%               Trial 3 in the second.
%               
%               If the field filter does not exist or is 
%               empty, all the data will be used.
%
%       DATA:   Specifies which statistical operation will 
%               be applied to which dimension of OUTPUT.
%               Like FILTER, it is a cell array with as many 
%               rows as the number of data subsets you want 
%               to compare.
%               Each row has the form
%                   {'Stat. op.' 'Specifier 1' ...
%                          'Stat. op.' 'Specifier 2' ...}
%               Statistical operation may be 'Mean',
%               'Median', 'Min', 'Max', 'Hist' or 'Pool'.
%               E.g., 
%
%                    plot_mode.data =
%                           {'Mean' 'Trial'; ...
%                           'Median' 'Trial'}
%
%               will calculate the average over all trials 
%               (which are left after filtering) for the 
%               first plot, and the median of all trials for
%               the second.
%               
%               If the field data does not exist or is 
%               empty, all the data will be pooled together
%               and the mean is calculated.
%
%       AUTOCOMBS: Instead of selecting each plot entry 
%               seperately by defining various rows in 
%               FILTER, AUTOCOMBS can be used to generate
%               automatically all the possible index 
%               combinations for a list of string 
%               specifiers. 
%               E.g., if #Groups = 2 and #Trials = 2 
%
%                    plot_mode.autocombs =
%                           {'Group Trial'}
%
%               will generate combinations 
%               (Group 1, Trial 1), (Group 1, Trial 2),  
%               (Group 2, Trial 1), (Group 2, Trial 2) 
%               and plot the corresponding data separately 
%               (in the same figure window).
%
%               Note that if FILTER is defined, combinations 
%               are only generated from the indices present
%               in FILTER. In this case, FILTER and DATA can
%               only have one row, and the statistical 
%               operation(s) in data will be applied to all 
%               combinations.   
%
%       LEGENDSTRING: Can be used to manually set the legend 
%               string. For each row in FILTER/each 
%               combination generated by autocombs, there
%               must be one row in LEGENDSTRING:
%                   plot_mode.legendstring =
%                           {'Mean over trials'; ...
%                           'Median of trials'}               
%
%       STATISTICS does not apply to all functions, since 
%               for many functions STATISTICS is fixed. It 
%               determines the statistical operation which 
%               is applied directly to the output of the 
%               core function. This only makes sense if 
%               the output is a numeric vector or a 
%               n -by- m cell of numeric vectors (e.g., a 1 
%               -by- m cell where each element of the cell  
%               corresponds to the bin of a histogram and 
%               contains all the values falling into this 
%               bin). 
%               The statistical operation can be 'Mean', 
%               'Median' and 'Hist' (the latter simply 
%               counts the number of values).
%               E.g., 
%
%                    plot_mode.statistics =
%                           'Mean';