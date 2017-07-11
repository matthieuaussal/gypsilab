function Mh = hmxChol(Mh)
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
%| Synopsis   : Cholesky factorization of H-Matrix                        |
%+========================================================================+

% H-Matrix (recursion)
if (Mh.typ == 0)
    % U11 -> M11
    Mh.chd{1} = hmxChol(Mh.chd{1});

    % U12 -> U11' \ M12
    Mh.chd{2} = hmxSolveLower(Mh.chd{1}',Mh.chd{2});
    
    % U21 -> 0
    Mh.chd{3}.dat = sparse(Mh.chd{3}.dim(1),Mh.chd{3}.dim(2));
    Mh.chd{3}.typ = 3;
    
    % M22 -> M22 - U12'*U12
    Mh.chd{4} = Mh.chd{4} - Mh.chd{2}' * Mh.chd{2};
    
    % U22 -> M22
    Mh.chd{4} = hmxChol(Mh.chd{4});
    
    % Fusion
    Mh = hmxFusion(Mh);

% Compressed leaf
elseif (Mh.typ == 1)
    error('hmxChol : unavailable case')
    
% Full leaf
elseif (Mh.typ == 2)
    Mh.dat = chol(Mh.dat);    

% Sparse leaf
elseif (Mh.typ == 3)
    Mh.dat = chol(Mh.dat);
   
% Unknown type
else
    error('hmxChol.m : unavailable case')
end
end
