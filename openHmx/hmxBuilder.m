function Mh = hmxBuilder(X,Y,green,tol)
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
%|    #    |   FILE       : hmxBuilder.m                                  |
%|    #    |   VERSION    : 0.31                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 25.11.2017                                    |
%| ( === ) |   SYNOPSIS   : Particles builder with low-rank approximation |
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
    Mh.typ = 0;    
end

% H-Matrix (recursion)
if (Mh.typ == 0)
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

    % H-Matrix (recursion)
    for i = 1:4
        Xi = X(Mh.row{i},:);
        Yi = Y(Mh.col{i},:);
        if isa(green,'function_handle')
            Mh.chd{i} = hmxBuilder(Xi,Yi,green,tol);
        elseif isnumeric(green)
            Mi = green(Mh.row{i},Mh.col{i});
            Mh.chd{i} = hmxBuilder(Xi,Yi,Mi,tol);
        else
            error('hmxBuilder.m : unavailable case')
        end
    end
    
    % Fusion
    Mh = hmxFusion(Mh);    

% Compressed leaf
elseif (Mh.typ == 1)
    % ACA with partial pivoting
    if isa(green,'function_handle')
        [A,B,flag] = hmxACA(X,Y,green,tol);

    % ACA with total pivoting
    elseif isnumeric(green) && ~issparse(green)
        [A,B,flag] = hmxACA(green,tol);
        
    % No compressor
    else
        flag = 0;
    end
    
    % Update
    if flag
        Mh.dat{1} = A;
        Mh.dat{2} = B;
    else
        Mh.typ = 2;
    end
end

% Full or sparse leaf
if (Mh.typ == 2)
    % Matrix
    if isa(green,'function_handle')
        [I,J]  = ndgrid(1:size(X,1),1:size(Y,1));
        Mh.dat = green(X(I,:),Y(J,:));
        Mh.dat = reshape(Mh.dat,size(X,1),size(Y,1));
    elseif isnumeric(green)
        Mh.dat = green;
    else
        error('hmxBuilder.m : unavailable case')
    end
    
    % Sparse matrix
    if issparse(Mh.dat)
        Mh.typ = 3;
    end
    
    % Full matrix recompression
    if (Mh.typ == 2)        
        [A,B,flag] = hmxSVD(Mh.dat,tol);
        if flag && (size(A,2) < 0.5*min(size(Mh.dat)))
            Mh.dat = {A,B};
            Mh.typ = 1;
        end
    end
end
end