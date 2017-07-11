function Mh = hmxBuilder(X,Y,green,tol)
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
%| Synopsis   : Particles builder with direct and low-rank separation     |
%+========================================================================+

% Initialisation
Mh = hmx(size(X,1),size(Y,1),tol);

% Rectangular box for X
Xmin = min(X,[],1);   
Xmax = max(X,[],1);
Xctr = 0.5*(Xmin+Xmax);
Xdgl = Xmax-Xmin;

% Rectangular box for Y
Ymin = min(Y,[],1);   
Ymax = max(Y,[],1);
Yctr = 0.5*(Ymin+Ymax);
Ydgl = Ymax-Ymin;

% Compression for separated boxes
if sum( abs(Yctr-Xctr) > max(Xdgl+Ydgl) ) 
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
        if isa(green,'function_handle')
            Mh.chd{i} = hmxBuilder(X(Mh.row{i},:),Y(Mh.col{i},:),green,tol);
        elseif isnumeric(green)
            Mh.chd{i} = hmxBuilder(X(Mh.row{i},:),Y(Mh.col{i},:), ...
                green(Mh.row{i},Mh.col{i}),tol);
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
        [A,B,flag] = hmxACA(X,Y,green,tol);
        
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
    
    % Full matrix recompression
    else
        rkMax      = ceil(sqrt(min(Mh.dim))-log(Mh.tol));
        [A,B,flag] = hmxRSVD(Mh.dat,Mh.tol,rkMax);
        if flag
            Mh.dat = {A,B};
            Mh.typ = 1;
        end
    end
end
end
