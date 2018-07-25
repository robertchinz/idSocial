function outvar=idSocial_auxiliaries_load(location)

% whos(var)
temp = load(location);
fn = fieldnames(temp);
var = hlp_deserialize(temp.(fn{1}));
outvar.(fn{1})=var;
