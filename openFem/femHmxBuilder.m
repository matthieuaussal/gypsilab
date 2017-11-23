function Mh = femHmxBuilder(Xdof,Ydof,Mx,X,green,Y,My,tol)
%+========================================================================+
%|                                                                        |
%|              OPENFEM - LIBRARY FOR FINITE ELEMENT METHOD               |
%|           openFem is part of the GYPSYLAB toolbox for Matlab           |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal & Francois Alouges (c) 2015-2017           |
%| PROPERTY  : Centre de Mathematiques Appliquees, Ecole polytechnique,   |
%| route de Saclay, 91128 Palaiseau, France. All rights reserved          |
%| LICENCE   :                                                            |
%| This program is a full open source software. It is distributed in the  |
%| hope that it will be useful, but without any warranty. This program is |
%| free only for public purpose: you can modify and/or adapt it only if   |
%| you share all your contributions and results, in open source and under |
%| the same licence. For other uses, you may contact us to obtain a       |
%| suitable license. Please acknowledge the gypsilab toolbox in programs  |
%| or publications in which you use it.                                   |
%| CONTACT   : matthieu.aussal@polytechnique.edu                          |
%|             francois.alouges@polytechnique.edu                         |
%| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
%|________________________________________________________________________|
%|   '&`   |                                                              |
%|    #    |   FILE       : femHmxBuilder.m                               |
%|    #    |   VERSION    : 0.31                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 25.11.2017                                    |
%| ( === ) |   SYNOPSIS   : Parallel particles builder with low-rank      |
%|  `---'  |                approximation and Finite Element Integration  |
%+========================================================================+

% Initialisation
Mh = hmx(size(Xdof,1),size(Ydof,1),tol);

% Centroids for each particles set
X0 = mean(Xdof,1);
Y0 = mean(Ydof,1);

% Unitary direction vector 
U = Y0 - X0;
if (norm(U) > 1e-6)
    U = U./norm(U);
else
    U = [1 0 0];
end

% Projection along the box separation
Xu = (Xdof-ones(size(Xdof,1),1)*X0) * U';
Yu = (Ydof-ones(size(Ydof,1),1)*X0) * U';

% Distances for particles X
Xr = sqrt(sum(  (Xdof-ones(size(Xdof,1),1)*X0).^2 , 2) );

% Compression for separated boxes
if (min(Yu) - max(Xu) > max(Xr))
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
        
        % Fem matrix subdivision and quadratures points indices
        [Mxchd,Ix] = femSubdivideCell(Mx,Ir,'left');
        [Mychd,Iy] = femSubdivideCell(My,Ic,'right');

        % Recursion
        Mh.chd{i} = femHmxBuilder(Xdof(Ir,:),Ydof(Ic,:),...
            Mxchd,X(Ix,:),green,Y(Iy,:),Mychd,tol);
    end
    
    % Fusion
    Mh = hmxFusion(Mh);    

    
% Compressed leaf
elseif (Mh.typ == 1)
    % ACA with partial pivoting
    if iscell(green)
        A = cell(1,3);
        B = cell(1,3);
        [A{1},B{1},flag1] = hmxACA(X,Y,green{1},tol);
        [A{2},B{2},flag2] = hmxACA(X,Y,green{2},tol);
        [A{3},B{3},flag3] = hmxACA(X,Y,green{3},tol);
        flag = flag1 * flag2 * flag3;
    else
        [A,B,flag] = hmxACA(X,Y,green,tol);
    end
        
    % Update
    if flag
        % Summation
        if iscell(Mx) && iscell(A) && iscell(My)
            A = [ Mx{1}*A{2} , - Mx{1}*A{3} , ...
                  Mx{2}*A{3} , - Mx{2}*A{1} , ...
                  Mx{3}*A{1} , - Mx{3}*A{2} ] ;
            B = [ B{2}*My{3} ; B{3}*My{2} ; ...
                  B{3}*My{1} ; B{1}*My{3} ; ...
                  B{1}*My{2} ; B{2}*My{1} ] ;

        elseif iscell(Mx) && iscell(A) && ~iscell(My)
            A = [Mx{1}*A{1} , Mx{2}*A{2} , Mx{3}*A{3}];
            B = [B{1}*My    ; B{2}*My    ; B{3}*My ];
            
        elseif iscell(Mx) && ~iscell(A) && iscell(My)
            A = [Mx{1}*A , Mx{2}*A , Mx{3}*A];
            B = [B*My{1} ; B*My{2} ; B*My{3} ];
            
        elseif ~iscell(Mx) && iscell(A) && iscell(My)
            A = [Mx*A{1}    , Mx*A{2}    , Mx*A{3}];
            B = [B{1}*My{1} ; B{2}*My{2} ; B{3}*My{3}];
        
        else
            A = Mx * A;
            B = B * My;
        end
        
        % Recompression
        [A,B]     = hmxQRSVD(A,B,tol);
        Mh.dat{1} = A;
        Mh.dat{2} = B;

    else
        Mh.typ = 2;
    end
end

% Full leaf
if (Mh.typ == 2)
    % Quadrature matrix
    Nx    = size(X,1);
    Ny    = size(Y,1);
    [I,J] = ndgrid(1:Nx,1:Ny);
    if iscell(green)
        Gxy{1} = reshape(green{1}(X(I(:),:),Y(J(:),:)),Nx,Ny);
        Gxy{2} = reshape(green{2}(X(I(:),:),Y(J(:),:)),Nx,Ny);
        Gxy{3} = reshape(green{3}(X(I(:),:),Y(J(:),:)),Nx,Ny);
    else
        Gxy = reshape(green(X(I(:),:),Y(J(:),:)),Nx,Ny);
    end
        
    % Matrix integration
    Mh.dat = femMultiplyCell(Mx,Gxy,My);
    
    % Recompression
    [A,B,flag] = hmxSVD(Mh.dat,tol);
    if flag && (size(A,2) < 0.5*min(size(Mh.dat)))
        Mh.dat = {A,B};
        Mh.typ = 1;
    end
end
end
