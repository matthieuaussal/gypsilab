function Bh = hmxSolveUpper(Bh,Uh)
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
%| Synopsis   : Solve upper H-Matrix system                               |
%+========================================================================+

%%% Security
if isstruct(Uh) && (Uh.typ == 0)
    if (Uh.chd{3}.typ ~= 3)
        error('hmxSolveUpper.m : left matrix is not an upper structure.')
    elseif ~isempty(find(Uh.chd{3}.dat,1))
        error('hmxSolveUpper.m : left matrix is not an upper structure.')
    end
end


%%% H-Matrix / H-Matrix -> H-Matrix
if isa(Bh,'hmx') && isa(Uh,'hmx')    
    % Check dimensions
    if (Bh.dim(2) ~= Uh.dim(2))
        error('hmxSolveUpper.m : matrix dimensions must agree.')
    end
    
    % H-Matrix / H-Matrix -> H-Matrix (recursion)
    if (Bh.typ == 0) && (Uh.typ == 0)   
        % X11 -> B11 / U11
        Bh.chd{1} = hmxSolveUpper(Bh.chd{1},Uh.chd{1});
        
        % X21 -> B21 / U11
        Bh.chd{3} = hmxSolveUpper(Bh.chd{3},Uh.chd{1});
        
        % X12 -> (B12 - X11*U12) / U22
        Bh.chd{2} = Bh.chd{2} - Bh.chd{1} * Uh.chd{2};
        Bh.chd{2} = hmxSolveUpper(Bh.chd{2},Uh.chd{4});
        
        % X22 -> (B22 - X21 * U12) / U22
        Bh.chd{4} = Bh.chd{4} - Bh.chd{3} * Uh.chd{2};
        Bh.chd{4} = hmxSolveUpper(Bh.chd{4},Uh.chd{4});
        
    % H-Matrix / Full -> Full
    elseif (Bh.typ == 0) && (Uh.typ == 2)
        B      = full(Bh) / Uh.dat;
        Bh     = hmx(size(B,1),size(B,2),Bh.tol);
        Bh.dat = B;
        Bh.typ = 2;
        
    % H-Matrix / Sparse -> Sparse
    elseif (Bh.typ == 0) && (Uh.typ == 3)
        B      = sparse(Bh) / Uh.dat;
        Bh     = hmx(size(B,1),size(B,2),Bh.tol);
        Bh.dat = B;
        Bh.typ = 3;

    % Compr / --- -> Compr
    elseif (Bh.typ == 1)
        Bh.dat = {Bh.dat{1} , hmxSolveUpper(Bh.dat{2},Uh)};
        
    % Full / --- -> Full
    elseif (Bh.typ == 2)
        Bh.dat = hmxSolveUpper(Bh.dat,Uh);
        
    % Sparse / --- -> Sparse
    elseif (Bh.typ == 3)
        Bh.dat = hmxSolveUpper(Bh.dat,Uh);

    else        
        error('hmxSolveUpper : unavailable case')
    end
        
    
%%% Matrix / H-Matrix -> Matrix   
elseif isa(Uh,'hmx')
    % Check dimensions
    if (size(Bh,2) ~= Uh.dim(2))
        error('hmxSolveUpper.m : matrix dimensions must agree.')
    end
    
    % H-Matrix (recursion)
    if (Uh.typ == 0)
        % X1 -> B1 / U11
        X1 = hmxSolveUpper(Bh(:,Uh.col{1}),Uh.chd{1});
        
        % X2 -> U22 / (B2 - X1*U12)
        X2 = Bh(:,Uh.col{2}) - X1*Uh.chd{2};
        X2 = hmxSolveUpper(X2,Uh.chd{4});
        
        % Final vector
        if issparse(X1)
            Bh = sparse(size(Bh,1),Uh.dim(2));
        else
            Bh = zeros(size(Bh,1),Uh.dim(2),class(X1));
        end
        Bh(:,Uh.row{1}) = X1;
        Bh(:,Uh.row{3}) = X2;
       
    % Compressed leaf
    elseif (Uh.typ == 1)
        error('hmxSolveUpper : unavailable case')

    % Full leaf
    elseif (Uh.typ == 2)
        Bh =  Bh / Uh.dat;
        
    % Sparse leaf
    elseif (Uh.typ == 3)
        Bh =  Bh / Uh.dat;
    
    % Unknown type    
    else
        error('hmxSolveUpper : unavailable case')
    end
  
    
    
%%% H-Matrix \ Matrix -> Matrix    
elseif isa(Bh,'hmx')
    % Check dimensions
    if (Bh.dim(1) ~= size(Uh,1))
        error('hmxSolveUpper.m : matrix dimensions must agree.')
    end
    
    % H-Matrix (recursion)
    if (Bh.typ == 0)
        % X2 -> B22 \ U2
        X2 = hmxSolveUpper(Bh.chd{4},Uh(Bh.row{3},:));
        
        % X1 -> B11 \ (U1 - B12*X2)
        X1 = Uh(Bh.row{1},:) - Bh.chd{2}*X2;
        X1 = hmxSolveUpper(Bh.chd{1},X1);

        % Uh = [X1 X2]
        if issparse(X1)
            Uh = sparse(Bh.dim(1),size(Uh,2));
        else
            Uh = zeros(Bh.dim(1),size(Uh,2),class(X1));
        end
        Uh(Bh.col{1},:) = X1;
        Uh(Bh.col{2},:) = X2;
        Bh              = Uh;
        
    % Compressed leaf
    elseif (Bh.typ == 1)
        error('hmxSolveUpper : unavailable case')

    % Full leaf
    elseif (Bh.typ == 2)
        Bh =  Bh.dat \ Uh;
    
    % Sparse leaf
    elseif (Bh.typ == 3)
        Bh =  Bh.dat \ Uh;
    
    % Unknown type    
    else
        error('hmxSolveUpper : unavailable case')
    end
  
    
    
%%% Unavailable
else
    error('hmxSolveUpper.m : unavailable case')
end
end
