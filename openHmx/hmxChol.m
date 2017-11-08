function Mh = hmxChol(Mh)
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
%|    #    |   FILE       : hmxChol.m                                     |
%|    #    |   VERSION    : 0.30                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 31.10.2017                                    |
%| ( === ) |   SYNOPSIS   : Cholesky factorization of H-Matrix            |
%|  `---'  |                                                              |
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
