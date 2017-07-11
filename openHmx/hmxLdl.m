function [Mh,Dh] = hmxLdl(Mh)
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
%| Synopsis   : LDLt factorization of H-Matrix                            |
%+========================================================================+

% H-Matrix (recursion)
if (Mh.typ == 0)
    % Digaonal initialisation
    Dh     = hmx(Mh.dim(1),Mh.dim(2),Mh.tol);
    Dh.row = Mh.row;
    Dh.col = Mh.col;
    Dh.typ = 0;
    
    % Nullify upper corner    
    Dh.chd{2}     = hmx(Mh.chd{2}.dim(1),Mh.chd{2}.dim(2),Mh.tol);
    Dh.chd{2}.dat = sparse(Dh.chd{2}.dim(1),Dh.chd{2}.dim(2));
    Dh.chd{2}.typ = 3;
    
    % Nullify lower corner    
    Dh.chd{3}     = hmx(Mh.chd{3}.dim(1),Mh.chd{3}.dim(2),Mh.tol);
    Dh.chd{3}.dat = sparse(Dh.chd{3}.dim(1),Dh.chd{3}.dim(2));
    Dh.chd{3}.typ = 3;
    
    % [L11,D11] -> M11
    [Mh.chd{1},Dh.chd{1}] = hmxLdl(Mh.chd{1});

    % L12 -> 0
    Mh.chd{2}.dat = sparse(Mh.chd{2}.dim(1),Mh.chd{2}.dim(2));
    Mh.chd{2}.typ = 3;
    
    % L21 -> M21 / U11
    Mh.chd{3} = hmxSolveUpper(Mh.chd{3},Dh.chd{1}*Mh.chd{1}');
    
    % M22 -> M22 - L21*U12
    Mh.chd{4} = Mh.chd{4} - Mh.chd{3} * Dh.chd{1} * Mh.chd{3}' ;
    
    % [L22,U22] -> M22
    [Mh.chd{4},Dh.chd{4}] = hmxLdl(Mh.chd{4});
    
    % Fusion
    Mh = hmxFusion(Mh);
    
% Compressed leaf
elseif (Mh.typ == 1)
    error('hmxLdl : unavailable case')
    
% Full leaf
elseif (Mh.typ == 2)
    Dh              = Mh;
    [Mh.dat,Dh.dat] = ldl(Mh.dat);      

% Sparse leaf
elseif (Mh.typ == 3)
    Dh      = Mh;
    [L,D,P] = ldl(Mh.dat);
    L       = P * L;
    Mh.dat  = L;
    Dh.dat  = D;    

% Unknown type
else
    error('hmxLdl.m : unavailable case')
end
end
