function idSocial_auxiliaries_save(location,varargin)


% Extract the variable values [From savefast RCH]
  vars = cell(size(varargin));
  for i = 1:numel(vars)
    vars{i} = evalin('caller', varargin{i});
  end 
  
%   for i = 1:numel(vars)
%       eval([varargin{i} ' = hlp_serialize(vars{i});']);
%   end
  for i = 1:numel(vars)
      eval([varargin{i} ' = vars{i};']);
  end

savefast(location,varargin{1});
