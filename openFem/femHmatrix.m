function Mh = femHmatrix(Xdof,Ydof,Mx,X,green,Y,My,tol)
%+========================================================================+
%|                                                                        |
%|                  OPENFEM, FINITE AND BOUNDARY ELEMENT                  |
%|              openFem is part of GYPSYLAB toolbox - v0.20               |
%|                                                                        |
%| Copyright (c) 20015-2017, Ecole polytechnique, all rights reserved.    |
%| Licence Creative Commons BY-NC-SA 4.0, Attribution, NonCommercial and  |
%| ShareAlike (see http://creativecommons.org/licenses/by-nc-sa/4.0/).    |
%| This software is the property from Centre de Mathematiques Appliquees  |
%| de l'Ecole polytechnique, route de Saclay, 91128 Palaiseau, France.    |
%|                                                            _   _   _   |
%| Please acknowledge the GYPSILAB toolbox in programs       | | | | | |  |
%| or publications in which you use the code. For openFem,    \ \| |/ /   |
%| we suggest as reference :                                   \ | | /    |
%| [1] : www.cmap.polytechnique.fr/~aussal/gypsilab             \   /     |
%|                                                               | |      |
%|_______________________________________________________________|_|______|
%| Author(s)  : Matthieu Aussal - CMAP, Ecole polytechnique               |
%| Creation   : 14.03.17                                                  |
%| Last modif : 21.06.17                                                  |
%| Synopsis   : Finite element builder with H-matrix                      |
%+========================================================================+

% Initialisation
Mh = hmx(size(Xdof,1),size(Ydof,1),tol);

% Rectangular box for Xdof
Xmin = min(Xdof,[],1);   
Xmax = max(Xdof,[],1);
Xctr = 0.5*(Xmin+Xmax);
Xdgl = Xmax-Xmin;

% Rectangular box for Ydof
Ymin = min(Ydof,[],1);   
Ymax = max(Ydof,[],1);
Yctr = 0.5*(Ymin+Ymax);
Ydgl = Ymax-Ymin;

% Compression for separated boxes
if sum( abs(Yctr-Xctr) > max(Xdgl+Ydgl) ) 
    Mh.typ = 1;

% Full computation for small box (stopping criterion)
elseif sum(Mh.dim < 100)
    Mh.typ = 2;

% H-Matrix (recursion)
else
    Mh.typ = 0;    
end

% H-Matrix (recursion)
if (Mh.typ == 0)
    % Subdivision for Xdof
    ind    = hmxSubdivide(Xdof);
    Mh.row = {find(ind),find(ind),find(~ind),find(~ind)};
    
    % Subdivision for Ydof
    ind    = hmxSubdivide(Ydof);
    Mh.col = {find(ind),find(~ind),find(ind),find(~ind)};
    
    % Single class
    if isa(Xdof,'single')
        for i = 1:4
            Mh.row{i} = single(Mh.row{i});
            Mh.col{i} = single(Mh.col{i});
        end
    end
    
    % H-Matrix (recursion)
    for i = 1:4
        % Dof indices
        Ir = Mh.row{i};
        Ic = Mh.col{i};
        
        % Quadrature indices
        if iscell(Mx) && iscell(My)
            Vx = 0;
            Vy = 0;
            for j = 1:3
                Vx = Vx + Mx{j}(:,Ir) * (1+rand(length(Ir),1));
                Vy = Vy + My{j}(:,Ic) * (1+rand(length(Ic),1));
            end
        else
           Vx = Mx(:,Ir) * (1+rand(length(Ir),1));
           Vy = My(:,Ic) * (1+rand(length(Ic),1));
        end
        Ix = find(Vx);
        Iy = find(Vy);
        
        % Security for empty matrix
        if (size(Ix,1) == 0)
            Ix = 1;
        end
        if (size(Iy,1) == 0)
            Iy = 1;
        end
        
        % Matrix subdivision
        if iscell(Mx) && iscell(My)
            Mxchd = cell(1,3);
            Mychd = cell(1,3);
            for j = 1:3
                Mxchd{j} = Mx{j}(Ix,Ir);
                Mychd{j} = My{j}(Iy,Ic);
            end
        else
            Mxchd = Mx(Ix,Ir);
            Mychd = My(Iy,Ic);
        end

        % Recursion
        Mh.chd{i} = femHmatrix(Xdof(Ir,:),Ydof(Ic,:),...
            Mxchd,X(Ix,:),green,Y(Iy,:),Mychd,tol);
    end
    
    % Fusion
    Mh = hmxFusion(Mh);    

    
% Compressed leaf
elseif (Mh.typ == 1)
    % ACA with partial pivoting
    [A,B,flag] = hmxACA(X,Y,green,tol);
        
    % Update
    if flag
        if iscell(Mx) && iscell(My)
            A     = [Mx{1}'*A , Mx{2}'*A , Mx{3}'*A];
            B     = [B*My{1}  ; B*My{2}  ; B*My{3} ];
            [A,B] = hmxQRSVD(A,B,tol);
        else
            A = Mx' * A;
            B = B * My;
        end
        Mh.dat{1} = A;
        Mh.dat{2} = B;
    else
        Mh.typ = 2;
    end
end


% Full leaf
if (Mh.typ == 2)
    % Quadrature matrix
    [I,J] = ndgrid(1:size(X,1),1:size(Y,1));
    Gxy   = green(X(I(:),:),Y(J(:),:));
    Gxy   = reshape(Gxy,size(X,1),size(Y,1));
    
    % Matrix integration
    if iscell(Mx) && iscell(My)
        Mh.dat = 0;
        for i = 1:3
            Mh.dat = Mh.dat + Mx{i}' * Gxy * My{i};
        end
    else
        Mh.dat = Mx' * Gxy * My;
    end
    
    % Full matrix recompression
    rkMax      = ceil(sqrt(min(Mh.dim))-log(Mh.tol));
    [A,B,flag] = hmxRSVD(Mh.dat,Mh.tol,rkMax);
    if flag
        Mh.dat = {A,B};
        Mh.typ = 1;
    end
end
end
