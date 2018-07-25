function sello=pos2sello(R)



x=-R:R;

X=repmat(x,[2*R+1 1]);
Y=X';

sello=sqrt(X.^2+Y.^2)<=R;



