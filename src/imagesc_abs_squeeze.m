% Shortcut for image dispaly
% 
% Author: F. Odille, IADI Nancy, 25/05/2011
% 

function imagesc_abs_squeeze(data, args)

if(nargin<2)
    imagesc( abs(squeeze(data)) )
else
    imagesc( abs(squeeze(data)), args )
end
