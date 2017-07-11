function Gxy = femGreenKernel(X,Y,green,k)
%+========================================================================+
%|                                                                        |
%|                  OPENFEM, FINITE AND BOUNDARY ELEMENT                  |
%|              openFem is part of GYPSYLAB toolbox - v0.20               |
%|                                                                        |
%| Copyright (c) 20015-2017, Ecole polytechnique, all rights reserved.    |
%| Licence Creative Commons BY-NC-SA 4.0, Attribution, NonCommercial and  |
%| ShareAlike (see http://creativecommons.org/licenses/by-nc-sa/4.0/).    |
%| This software is the property from Centre de Mathematiques Appliquees  |
%| de l'Ecole polytechnique, route de Saclay, 91128 Palaiseau, France.    |
%|                                                            _   _   _   |
%| Please acknowledge the GYPSILAB toolbox in programs       | | | | | |  |
%| or publications in which you use the code. For openFem,    \ \| |/ /   |
%| we suggest as reference :                                   \ | | /    |
%| [1] : www.cmap.polytechnique.fr/~aussal/gypsilab             \   /     |
%|                                                               | |      |
%|_______________________________________________________________|_|______|
%| Author(s)  : Matthieu Aussal - CMAP, Ecole polytechnique               |
%| Creation   : 14.03.17                                                  |
%| Last modif : 21.06.17                                                  |
%| Synopsis   : Classical green kernel function                           |
%+========================================================================+

% Security
if (size(X,2) ~= 3) || (size(Y,2) ~= 3)
    error('femGreenKernel.m : unavailable case')
end
if isempty(k)
    k = 0;
end

% Distances between particles
Rxy = sqrt( ...
    (X(:,1) - Y(:,1)).^2 + ...
    (X(:,2) - Y(:,2)).^2 + ...
    (X(:,3) - Y(:,3)).^2 );

% For empty wave-number
if isempty(k)
    k = 0;
end

% Green kernel definition
if strcmp(green,'[1/r]')
    Gxy = 1./Rxy;   
    
elseif strcmp(green,'[exp(ikr)/r]')
    Gxy = exp(1i*k*Rxy)./Rxy;          

elseif strcmp(green(1:end-1),'gradx[1/r]')    
    j = str2double(green(end));
    Gxy = - (X(:,j)-Y(:,j)) ./ (Rxy.^3);
    
elseif strcmp(green(1:end-1),'grady[1/r]')     
    j = str2double(green(end));
    Gxy = (X(:,j)-Y(:,j)) ./ (Rxy.^3);    
    
elseif strcmp(green(1:end-1),'gradx[exp(ikr)/r]')  
    j = str2double(green(end));
    Gxy = (1i*k - 1./Rxy) .* exp(1i*k.*Rxy) .* ...
        (X(:,j)-Y(:,j)) ./ (Rxy.^2);
    
elseif strcmp(green(1:end-1),'grady[exp(ikr)/r]')    
    j = str2double(green(end));
    Gxy = - (1i*k - 1./Rxy) .* exp(1i*k.*Rxy) .* ...
        (X(:,j)-Y(:,j)) ./ (Rxy.^2);
    
elseif strcmp(green(1:end-2),'[ij/r+rirj/r^3]')        
    i = str2double(green(end-1));
    j = str2double(green(end));
    Gxy = (i==j)./Rxy + (X(:,i)-Y(:,i)).*(X(:,j)-Y(:,j))./(Rxy.^3);
    
elseif strcmp(green(1:end-3),'[rirjrk/r^5]')      
    i = str2double(green(end-2));
    j = str2double(green(end-1));
    k = str2double(green(end));    
    Gxy = (X(:,i)-Y(:,i)).*(X(:,j)-Y(:,j)).*(X(:,k)-Y(:,k))./(Rxy.^5);
    
else
    error('Error in hmxGreenKernel.m : unknown green kernel')
end

% Singularity
if strcmp(green,'[exp(ikr)/r]')
    Gxy(Rxy<1e-6) = 0 + 1i*k;
else
    Gxy(Rxy<1e-6) = 0;
end

end
