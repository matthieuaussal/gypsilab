function Mh = hmxTree(X,Y,tol)
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
%|    #    |   FILE       : hmxTree.m                                     |
%|    #    |   VERSION    : 0.30                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 31.10.2017                                    |
%| ( === ) |   SYNOPSIS   : Binary tree with median repartition           |
%|  `---'  |                                                              |
%+========================================================================+

% Initialisation
Mh = hmx(size(X,1),size(Y,1),tol);

% Centroids for each particles set
X0 = mean(X,1);
Y0 = mean(Y,1);

% Unitary direction vector 
U = Y0 - X0;
if (norm(U) > 1e-6)
    U = U./norm(U);
else
    U = [1 0 0];
end

% Projection along the box separation
Xu = (X-ones(size(X,1),1)*X0) * U';
Yu = (Y-ones(size(Y,1),1)*X0) * U';

% Distances for particles X
Xr = sqrt(sum(  (X-ones(size(X,1),1)*X0).^2 , 2) );

% Compression for separated boxes
if (min(Yu) - max(Xu) > max(Xr))
    Mh.typ = 1;

% Full computation for small box (stopping criterion)
elseif sum(Mh.dim < 100)
    Mh.typ = 2;

% H-Matrix
else
    % Subdivision for X
    ind    = hmxSubdivide(X);
    Mh.row = {find(ind),find(ind),find(~ind),find(~ind)};
    
    % Subdivision for Y
    ind    = hmxSubdivide(Y);
    Mh.col = {find(ind),find(~ind),find(ind),find(~ind)};
    
    % Single class
    if isa(X,'single')
        for i = 1:4
            Mh.row{i} = single(Mh.row{i});
            Mh.col{i} = single(Mh.col{i});
        end
    end

    % Leaves recursion
    typ = zeros(1,4);
    for i = 1:4
        Mh.chd{i} = hmxTree(X(Mh.row{i},:),Y(Mh.col{i},:),tol);
        typ(i)    = Mh.chd{i}.typ;
    end
    
    % Low-rank fusion
    if (sum(typ==1) == 4)
        Mh     = hmx(size(X,1),size(Y,1),tol);
        Mh.typ = 1;

    % H-Matrix    
    else    
        Mh.typ = 0;
    end
end
end
