function combos = idSocial_auxiliaries_allcombs(vec1, vec2, vec3, vec4, vec5, vec6)
combos = [];
if nargin==2
[p,q] = meshgrid(vec1, vec2);
combos = [p(:) q(:)];
elseif nargin==3

    [p,q,r] = meshgrid(vec1, vec2, vec3);
    combos = [p(:) q(:) r(:)];
elseif nargin==6

    [p,q,r,s,t,u] = ndgrid(vec1, vec2, vec3, vec4, vec5, vec6);
    combos = [p(:) q(:) r(:) s(:) t(:) u(:)];
end