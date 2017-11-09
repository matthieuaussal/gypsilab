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
%|    #    |   VERSION    : 0.30                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 31.10.2017                                    |
%| ( === ) |   SYNOPSIS   : Parallel particles builder with low-rank      |
%|  `---'  |                approximation                                 |
%+========================================================================+

% Recursive binary tree
Mh = hmxTree(X,Y,tol);

% Recursive leaves blocks
Ml = hmxLeafOut(Mh);

% Randomize block repartition
Iprm = randperm(length(Ml));
Ml   = Ml(Iprm);

% Parallel computing for partial pivoting only
if isnumeric(green)
    Nlab = 1;
else
    Nlab = length(Composite);
end

% Parallel computation for partial pivoting
parfor (n = 1:length(Ml),Nlab)
    % Leaf data
    I   = Ml{n}{1};
    J   = Ml{n}{2};
    dat = Ml{n}{3};
    typ = Ml{n}{4};
        
    % ACA for compressed leaf
    if (typ == 1)
        % Partial pivoting
        if isa(green,'function_handle')
            [A,B,flag] = hmxACA(X(I,:),Y(J,:),green,tol);
            
        % Total pivoting
        elseif isnumeric(green) && ~issparse(green)
            [A,B,flag] = hmxACA(green(I,J),tol);

        % No compression
        else
            A    = [];
            B    = [];
            flag = 0;
        end
        
        % Update
        if flag
            dat = {A,B};
        else
            typ = 2;
        end
    end
    
    % No compression for full or sparse leaf 
    if (typ == 2)
        if isa(green,'function_handle')
            [idx,jdx] = ndgrid(I,J);
            dat       = green(X(idx,:),Y(jdx,:));
            dat       = reshape(dat,length(I),length(J));
        elseif isnumeric(green)
            dat = green(I,J);
        else
            error('hmxBuilder.m : unavailable case')
        end
    end
    
    % Sparse matrix
    if issparse(dat)
        typ = 3;
    end
    
    % Full matrix recompression
    if (typ == 2)        
        [A,B,flag] = hmxSVD(dat,tol);
        if flag && (size(A,2) < 0.5*min(size(dat)))
            dat = {A,B}
            typ = 1;
        end
    end
    
    % Update
    Ml{n}{3} = dat;
    Ml{n}{4} = typ;
    Ml{n}{5} = tol;
end

% Reorder block
Ml(Iprm) = Ml;

% Block repartition
Mh = hmxLeafIn(Mh,Ml);
end
