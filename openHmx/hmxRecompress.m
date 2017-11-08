function Mh = hmxRecompress(Mh,tol)
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
%|    #    |   FILE       : hmxRecompress.m                               |
%|    #    |   VERSION    : 0.30                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 31.10.2017                                    |
%| ( === ) |   SYNOPSIS   : Recompreesion of H-Matrix with new accuracy   |
%|  `---'  |                                                              |
%+========================================================================+

% H-Matrix (recursion)
if (Mh.typ == 0)
    for i = 1:4
        Mh.chd{i} = hmxRecompress(Mh.chd{i},tol);
        Mh.tol    = max(tol,Mh.tol);
    end
    Mh = hmxFusion(Mh);
    
% Compressed leaf
elseif (Mh.typ == 1)
    if (Mh.tol < tol)
        [A,B]     = hmxQRSVD(Mh.dat{1},Mh.dat{2},tol);
        Mh.dat{1} = A;
        Mh.dat{2} = B;
        Mh.tol    = tol;
    end
    
% Others leaves
else
    Mh.tol = max(tol,Mh.tol); 
end
end
