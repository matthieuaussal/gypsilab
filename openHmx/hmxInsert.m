function Mh = hmxInsert(Mh,S)
%+========================================================================+
%|                                                                        |
%|         OPENHMX - LIBRARY FOR H-MATRIX COMPRESSION AND ALGEBRA         |
%|           openHmx is part of the GYPSILAB toolbox for Matlab           |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal (c) 2017-2018.                             |
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
%|    #    |   FILE       : hmxInsert.m                                   |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Insert sparse matrix inside H-Matrix          |
%|  `---'  |                                                              |
%+========================================================================+

% Check input
if ~issparse(S)
    error('hmxInsert.m : unavailable case')
end

% H-Matrix (recursion)
if (Mh.typ == 0)
    for i = 1:4
        Mh.chd{i} = hmxInsert(Mh.chd{i},S(Mh.row{i},Mh.col{i}));
    end
    
% Compressed leaf
elseif (Mh.typ == 1) 
    if (nnz(S) > 0)
        error('hmxInsert.m : unavailable case')
    end
    
% Full leaf
elseif (Mh.typ == 2)
    [I,J,V] = find(S);
    Mh.dat(sub2ind(size(Mh.dat),I,J)) = V;
    
% Unknown type
else
    error('hmxInsert.m : unavailable case')
end
end
