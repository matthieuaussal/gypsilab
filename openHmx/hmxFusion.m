function Mh = hmxFusion(Mh)
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
%|    #    |   FILE       : hmxFusion.m                                   |
%|    #    |   VERSION    : 0.30                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 31.10.2017                                    |
%| ( === ) |   SYNOPSIS   : Fusion and recompress for full, sparse and    |
%|  `---'  |                low-rank leaves                               |
%+========================================================================+

% H-Matrix
if (Mh.typ == 0)
    % Leaves type
    typ = zeros(1,4);
    for i = 1:4
        typ(i) = Mh.chd{i}.typ;
    end
    
    % Low-rank fusion
    if (sum(typ==1) == 4)
        % Rank for each leaf
        rk = zeros(1,4);
        nk = 0;
        for i = 1:4
            rk(i) = size(Mh.chd{i}.dat{1},2);
            nk    = nk + sum(Mh.chd{i}.dim)*rk(i); 
        end
        
        % Low rank matrix
        A = zeros(Mh.dim(1),sum(rk),class(Mh.chd{1}.dat{1}));
        B = zeros(sum(rk),Mh.dim(2),class(Mh.chd{1}.dat{2}));
        j = 0;
        for i = 1:4
            A(Mh.row{i},j+(1:rk(i))) = Mh.chd{i}.dat{1};
            B(j+(1:rk(i)),Mh.col{i}) = Mh.chd{i}.dat{2};            
            j = j + rk(i);
        end
        [A,B] = hmxQRSVD(A,B,Mh.tol);
        
        % Update
        if (sum(Mh.dim)*size(A,2) < nk)
            Mh     = hmx(Mh.dim(1),Mh.dim(2),Mh.tol);
            Mh.dat = {A,B};
            Mh.typ = 1;
        end
    end
    
    % Sparse fusion
    if (sum(typ==3) == 4)
        % Emptyness for each leaf
        emp = 0;
        for i = 1:4
            emp = emp + isempty(find(Mh.chd{i}.dat,1));
        end
        
        % Update empty
        if (emp == 4)
            Mh     = hmx(Mh.dim(1),Mh.dim(2),Mh.tol);
            Mh.dat = sparse(Mh.dim(1),Mh.dim(2));
            Mh.typ = 3;
        end
        
        % Update sparse
        if (emp == 0)
            M      = sparse(Mh);
            Mh     = hmx(Mh.dim(1),Mh.dim(2),Mh.tol);
            Mh.dat = M;
            Mh.typ = 3;
        end
    end
    
% Unvalid case 
else
    error('hmxFusion.m : unavailable case')
end
end
