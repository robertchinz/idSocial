function cmInteract = idSocial_auxiliaries_generateInteractionColormap(clims,size_cm)

if nargin<1 || isempty(clims)
    clims = get(gca,'CLim');
end

if nargin<2 || isempty(size_cm)
    size_cm = 256;
end

if clims(1)<0 && clims(2)>0
    no_red = (size_cm * clims(2))/(clims(2)-clims(1));
    no_blue = (size_cm * abs(clims(1)))/(clims(2)-clims(1));
    
%     cmInteractRed = colormap4InteractionMaps(129:256,:);
%     cmInteractBlue = colormap4InteractionMaps(1:128,:);
    
    cmInteractRed = [1.0000, 1.0000, 0; ...
        0.5000, 0, 0];
    
    cmInteractBlue = [0, 0, 0.5156; ...
        0.0156, 1.0000, 0.9844];
    
    cmInteractRed2 = interp1([1,no_red],cmInteractRed,1:no_red);
    cmInteractBlue2 = interp1([1,no_blue],cmInteractBlue,1:no_blue);
    
    cmInteract = vertcat(cmInteractBlue2,cmInteractRed2);
else
    warning([mfilename ': Color limits must enclose 0.'])
end