function Ml = hmxLeafOut(varargin)
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
%|    #    |   FILE       : hmxLeafOut.m                                  |
%|    #    |   VERSION    : 0.30                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 31.10.2017                                    |
%| ( === ) |   SYNOPSIS   : Convert H-Matrix format to leaf format        |
%|  `---'  |                                                              |
%+========================================================================+

% Input analysis
Mh = varargin{1};
if (nargin == 1)
    I = (1:Mh.dim(1))';
    J = (1:Mh.dim(2))';
elseif (nargin == 3)
    I = varargin{2};
    J = varargin{3};
else
    error('hmxLeafOut.m : unavailable case.')
end

% H-Matrix (recursion)
if (Mh.typ == 0)
    Ml = cell(0,1);
    for i = 1:4
        tmp       = hmxLeafOut(Mh.chd{i},I(Mh.row{i}),J(Mh.col{i}));
        ind       = size(Ml,1) + (1:length(tmp));
        Ml(ind,1) = tmp(:);
    end
    
% Leaf
else
    Ml{1} = {I,J,Mh.dat,Mh.typ,Mh.tol};
end
end
