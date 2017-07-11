function Mh = hmxInv(Mh)
%+========================================================================+
%|                                                                        |
%|               OPENHMX, H-MATRIX COMPRESSION AND ALGEBRA                |
%|              openHmx is part of GYPSYLAB toolbox - v0.20               |
%|                                                                        |
%| Copyright (c) 20015-2017, Ecole polytechnique, all rights reserved.    |
%| Licence Creative Commons BY-NC-SA 4.0, Attribution, NonCommercial and  |
%| ShareAlike (see http://creativecommons.org/licenses/by-nc-sa/4.0/).    |
%| This software is the property from Centre de Mathematiques Appliquees  |
%| de l'Ecole polytechnique, route de Saclay, 91128 Palaiseau, France.    |
%|                                                            _   _   _   |
%| Please acknowledge the GYPSILAB toolbox in programs       | | | | | |  |
%| or publications in which you use the code. For openHmx,    \ \| |/ /   |
%| we suggest as reference :                                   \ | | /    |
%| [1] : www.cmap.polytechnique.fr/~aussal/gypsilab             \   /     |
%| [2] : 13th International Conference on Mathematical           | |      |
%| and Numerical Aspects of Wave Propagation, University of      | |      |
%| Minnesota, may 2017. "OpenHmX, an open-source H-Matrix        | |      |
%| toolbox in Matlab".                                           | |      |
%|_______________________________________________________________|_|______|
%| Author(s)  : Matthieu Aussal - CMAP, Ecole polytechnique               |
%| Creation   : 14.03.17                                                  |
%| Last modif : 21.06.17                                                  |
%| Synopsis   : Inversion of H-Matrix based on Schur complement           |
%+========================================================================+

% H-Matrix (bloc recursion)
if (Mh.typ == 0)
    % Am1 -> M11 
    Mh.chd{1} = hmxInv(Mh.chd{1});
    
    % - Am1 * B -> X12
    X12 = - Mh.chd{1} * Mh.chd{2};
    
    % C * Am1 -> X21
    X21 = Mh.chd{3} * Mh.chd{1};
    
    % D - C * Am1 * B = S -> M22  (schur complement)
    Mh.chd{4} = Mh.chd{4} + Mh.chd{3} * X12;
    
    % S^-1 -> M22
    Mh.chd{4} = hmxInv(Mh.chd{4});
    
    % - Am1 * B * S -> M12
    Mh.chd{2} = X12 * Mh.chd{4};
    
    % Am1 + Am1 * B * S * C * Am1 -> M11
    Mh.chd{1} = Mh.chd{1} - Mh.chd{2} * X21;
    
    % - S * C * Am1 -> M21
    Mh.chd{3} = - Mh.chd{4} * X21;
    
    % Fusion
    Mh = hmxFusion(Mh);
    
% Compressed leaf    
elseif (Mh.typ == 1)
    error('hmxInv : unavailable case')
    
% Full leaf    
elseif (Mh.typ == 2)
    Mh.dat = inv(Mh.dat);

% Sparse leaf    
elseif (Mh.typ == 3)
    Mh.dat = inv(Mh.dat);

% Unknown type
else
    error('hmxInv.m : unavailable case')
end
end
