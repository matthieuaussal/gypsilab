function Mh = hmxBuilderFem(Xdof,Ydof,Mx,X,green,Y,My,tol)
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
%|    #    |   FILE       : hmxBuilderFem.m                               |
%|    #    |   VERSION    : 0.30                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 31.10.2017                                    |
%| ( === ) |   SYNOPSIS   : Parallel particles builder with low-rank      |
%|  `---'  |                approximation and Finite Element Integration  |
%+========================================================================+

% Recursive binary tree
Mh = hmxTree(Xdof,Ydof,tol);

% Recursive leaves blocks
Mleaf = hmxLeafOut(Mh);

% Randomize block repartition
Iprm  = randperm(length(Mleaf));
Mleaf = Mleaf(Iprm);

% Parallel computation with partial pivoting
parfor n = 1:length(Mleaf)
    % Leaf data
    I   = Mleaf{n}{1};
    J   = Mleaf{n}{2};
    dat = Mleaf{n}{3};
    typ = Mleaf{n}{4};
    
    % Left integration indices
    if iscell(Mx)
        Vx = (1+rand(1,length(I))) * Mx{1}(I,:) + ...
            (1+rand(1,length(I))) * Mx{2}(I,:) + ...
            (1+rand(1,length(I))) * Mx{3}(I,:);
    else
        Vx = (1+rand(1,length(I))) * Mx(I,:);
    end
    Ix = find(Vx);
    if isempty(Ix)
        Ix = [1;2]';
    end
    
    % Left integration matrix
    if iscell(Mx)
        Ml{1} = Mx{1}(I,Ix);
        Ml{2} = Mx{2}(I,Ix);
        Ml{3} = Mx{3}(I,Ix);
    else
        Ml = Mx(I,Ix);
    end
    
    % Right integration indices
    if iscell(My)
        Vy = My{1}(:,J) * (1+rand(length(J),1)) + ...
            My{2}(:,J) * (1+rand(length(J),1)) + ...
            My{3}(:,J) * (1+rand(length(J),1)) ;
    else
        Vy = My(:,J) * (1+rand(length(J),1));
    end
    Iy = find(Vy);
    if isempty(Iy)
        Iy = [1;2]';
    end
    
    % Right integration matrix
    if iscell(My)
        Mr{1} = My{1}(Iy,J);
        Mr{2} = My{2}(Iy,J);
        Mr{3} = My{3}(Iy,J);
    else
        Mr = My(Iy,J);
    end

    % Compressed leaf
    if (typ == 1)
        % ACA with partial pivoting
        if iscell(green)
            A = cell(1,3);
            B = cell(1,3);
            [A{1},B{1},flag1] = hmxACA(X(Ix,:),Y(Iy,:),green{1},tol);
            [A{2},B{2},flag2] = hmxACA(X(Ix,:),Y(Iy,:),green{2},tol);
            [A{3},B{3},flag3] = hmxACA(X(Ix,:),Y(Iy,:),green{3},tol);
            flag = flag1 * flag2 * flag3;
        else
            [A,B,flag] = hmxACA(X(Ix,:),Y(Iy,:),green,tol);
        end
        
        % Compressed leaf
        if flag
            % Finite element integration
            if iscell(Ml) && iscell(A) && iscell(Mr)
                A = [ Ml{1}*A{2} , - Ml{1}*A{3} , ...
                    Ml{2}*A{3} , - Ml{2}*A{1} , ...
                    Ml{3}*A{1} , - Ml{3}*A{2} ] ;
                B = [ B{2}*Mr{3} ; B{3}*Mr{2} ; ...
                    B{3}*Mr{1} ; B{1}*Mr{3} ; ...
                    B{1}*Mr{2} ; B{2}*Mr{1} ] ;
                
            elseif iscell(Ml) && iscell(A) && ~iscell(Mr)
                A = [Ml{1}*A{1} , Ml{2}*A{2} , Ml{3}*A{3}];
                B = [B{1}*Mr    ; B{2}*Mr    ; B{3}*Mr ];
                
            elseif iscell(Ml) && ~iscell(A) && iscell(Mr)
                A = [Ml{1}*A , Ml{2}*A , Ml{3}*A];
                B = [B*Mr{1} ; B*Mr{2} ; B*Mr{3} ];
                
            elseif ~iscell(Ml) && iscell(A) && iscell(Mr)
                A = [Ml*A{1}    , Ml*A{2}    , Ml*A{3}];
                B = [B{1}*Mr{1} ; B{2}*Mr{2} ; B{3}*Mr{3}];
                
            else
                A = Ml * A;
                B = B * Mr;
            end
            
            % Recompression
            [A,B] = hmxQRSVD(A,B,tol);
            
            % Update
            dat = {A,B};

        else
            typ = 2;
        end
    end
    
    % Full leaf
    if (typ == 2)
        % Full gridding on quadrature points
        [idx,jdx] = ndgrid(Ix,Iy);
        
        % Quadrature matrix
        if iscell(green)
            Gxy{1} = reshape(green{1}(X(idx(:),:),Y(jdx(:),:)),length(Ix),length(Iy));
            Gxy{2} = reshape(green{2}(X(idx(:),:),Y(jdx(:),:)),length(Ix),length(Iy));
            Gxy{3} = reshape(green{3}(X(idx(:),:),Y(jdx(:),:)),length(Ix),length(Iy));
        else
            Gxy = reshape(green(X(idx(:),:),Y(jdx(:),:)),length(Ix),length(Iy));
        end
        
        % Matrix integration
        dat = femMultiplyCell(Ml,Gxy,Mr);
        
        % Recompression        
        [A,B,flag] = hmxSVD(dat,tol);
        if flag && (size(A,2) < 0.5*min(size(dat)))
            dat = {A,B}
            typ = 1;
        end
    end
    
    % Update
    Mleaf{n}{3} = dat;
    Mleaf{n}{4} = typ;
    Mleaf{n}{5} = tol;
end

% Reorder block
Mleaf(Iprm) = Mleaf;

% Block repartition
Mh = hmxLeafIn(Mh,Mleaf);
end
