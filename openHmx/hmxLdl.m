function [Mh,Dh] = hmxLdl(Mh)
%+========================================================================+
%|                                                                        |
%|         OPENHMX - LIBRARY FOR H-MATRIX COMPRESSION AND ALGEBRA         |
%|           openHmx is part of the GYPSILAB toolbox for Matlab           |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal (c) 2015-2017.                             |
%| PROPERTY  : Centre de Mathematiques Appliquees, Ecole polytechnique,   |
%| route de Saclay, 91128 Palaiseau, France. All rights reserved.         |
%| LICENCE   : This program is free software, distributed in the hope that|
%| it will be useful, but WITHOUT ANY WARRANTY. Natively, you can use,    |
%| redistribute and/or modify it under the terms of the GNU General Public|
%| License, as published by the Free Software Foundation (version 3 or    |
%| later,  http://www.gnu.org/licenses). For private use, dual licencing  |
%| is available, please contact us to activate a "pay for remove" option. |
%| CONTACT   : matthieu.aussal@polytechnique.edu                          |
%| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
%|                                                                        |
%| Please acknowledge the gypsilab toolbox in programs or publications in |
%| which you use it.                                                      |
%|________________________________________________________________________|
%|   '&`   |                                                              |
%|    #    |   FILE       : hmxLdl.m                                      |
%|    #    |   VERSION    : 0.30                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 31.10.2017                                    |
%| ( === ) |   SYNOPSIS   : LDLt factorization of H-Matrix                |
%|  `---'  |                                                              |
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
