% 10 feb 10 Corrección: Evito que la media se haga sólo por un lado en los
% extremos. En vez de eso, voy reduciendo igualmente por ambos lados
% (cuando n_pasos sea par, siempre habrá un paso más por detrás que por
% delante)
% APE 27 abr 10
% Modified using meanangle.m from Matlab Central FEX RCH, 9/3/2012
function vector=mediamovil4angles(vector,n_pasos)
disp([ '[' mfilename '] ' 'Caution: This function is not working properly. Check why!'])
vector=vector(:)';
l=length(vector);
despl=-floor((n_pasos-1)/2):ceil((n_pasos-1)/2);
matriz=NaN(n_pasos,l);
for c_despl=1:length(despl)
    matriz(c_despl,max([1 despl(c_despl)+1]):min([l l+despl(c_despl)]))=vector(max([1 -despl(c_despl)+1]):min([l l-despl(c_despl)]));
end % c_despl

% Simetrizamos los extemos para que no haya promedios raros en los extremos
for c=1:ceil(length(despl)/2)-1
        matriz(c,isnan(matriz(end-c+1,:)))=NaN;
        matriz(end-c+1,isnan(matriz(c,:)))=NaN;
end % c
vector=meanangle(matriz);
% vector=mean(matriz);