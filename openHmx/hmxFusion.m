function Mh = hmxFusion(Mh)
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
%| Synopsis   : Fusion and recompress for full, sparse and low-ran  k     |
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
        for i = 1:4
            rk(i) = size(Mh.chd{i}.dat{1},2);
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
        Mh     = hmx(Mh.dim(1),Mh.dim(2),Mh.tol);
        Mh.dat = {A,B};
        Mh.typ = 1;
    end
       
    % Full fusion
    if (sum(typ==2) == 4)
        M      = full(Mh);
        Mh     = hmx(Mh.dim(1),Mh.dim(2),Mh.tol);
        Mh.dat = M;
        Mh.typ = 2;
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
    error('hmxBuilder.m : unavailable case')
end
end
