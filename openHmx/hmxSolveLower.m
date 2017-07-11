function Bh = hmxSolveLower(Lh,Bh)
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
%| Synopsis   : Solve lower H-Matrix system                               |
%+========================================================================+

%%% Security
if isstruct(Lh) && (Lh.typ == 0)
    if (Lh.chd{2}.typ ~= 3)
        error('hmxSolveLower.m : left matrix is not a lower structure.')
    elseif ~isempty(find(Lh.chd{2}.dat,1))
        error('hmxSolveLower.m : left matrix is not a lower structure.')
    end
end


%%% H-Matrix \ H-Matrix -> H-Matrix
if isa(Lh,'hmx') && isa(Bh,'hmx')    
    % Check dimensions
    if (Lh.dim(1) ~= Bh.dim(1))
        error('hmxSolveLower.m : matrix dimensions must agree.')
    end
    
    % H-Matrix \ H-Matrix -> H-Matrix (recursion)
    if (Lh.typ == 0) && (Bh.typ == 0) 
        % X11 -> L11 \ B11
        Bh.chd{1} = hmxSolveLower(Lh.chd{1},Bh.chd{1});
        
        % X12 -> L11 \ B12
        Bh.chd{2} = hmxSolveLower(Lh.chd{1},Bh.chd{2});
        
        % X21 -> L22 \ (B21 - L21*X11)
        Bh.chd{3} = Bh.chd{3} - Lh.chd{3} * Bh.chd{1};
        Bh.chd{3} = hmxSolveLower(Lh.chd{4},Bh.chd{3});
        
        % X22 -> L22 \ (B22 - L21 * X12)
        Bh.chd{4} = Bh.chd{4} - Lh.chd{3} * Bh.chd{2};
        Bh.chd{4} = hmxSolveLower(Lh.chd{4},Bh.chd{4});
        
    % Full \ H-Matrix -> Full
    elseif (Lh.typ == 2) && (Bh.typ == 0)
        B      = Lh.dat \ full(Bh);
        Bh     = hmx(size(B,1),size(B,2),Bh.tol);
        Bh.dat = B;
        Bh.typ = 2;
    
    % Sparse \ H-Matrix -> Sparse
    elseif (Lh.typ == 3) && (Bh.typ == 0)
        B      = Lh.dat \ sparse(Bh);
        Bh     = hmx(size(B,1),size(B,2),Bh.tol);
        Bh.dat = B;
        Bh.typ = 3;        
        
    % --- \ Compr -> Compr
    elseif (Bh.typ == 1)
        Bh.dat = { hmxSolveLower(Lh,Bh.dat{1}) , Bh.dat{2} };
                
    % --- \ Full -> Full
    elseif (Bh.typ == 2)
        Bh.dat = hmxSolveLower(Lh,Bh.dat);
    
    % --- \ Sparse -> Sparse
    elseif (Bh.typ == 3)
        Bh.dat = hmxSolveLower(Lh,Bh.dat);
        
    else
        error('hmxSolveLower : unavailable case')
    end
    
    
%%% H-Matrix \ Matrix -> Matrix   
elseif isa(Lh,'hmx')
    % Check dimensions
    if (Lh.dim(1) ~= size(Bh,1))
        error('hmxhmxSolveLower.m : matrix dimensions must agree.')
    end
    
    % H-Matrix (recursion)
    if (Lh.typ == 0)
        % X1 -> L11 \ B1
        X1 = hmxSolveLower(Lh.chd{1},Bh(Lh.row{1},:));
        
        % X2 -> L22 \ (B2 - L21*X1)
        X2 = Bh(Lh.row{3},:) - Lh.chd{3}*X1;
        X2 = hmxSolveLower(Lh.chd{4},X2);
        
        % Bh = [X1 X2]
        if issparse(X1)
            Bh = sparse(Lh.dim(1),size(Bh,2));
        else
            Bh = zeros(Lh.dim(1),size(Bh,2),class(X1));
        end
        Bh(Lh.col{1},:) = X1;
        Bh(Lh.col{2},:) = X2;
              
    % Compressed leaf
    elseif (Lh.typ == 1)
        error('hmxSolveLower : unavailable case')
        
    % Full leaf
    elseif (Lh.typ == 2)
        Bh = Lh.dat \ Bh;

    % Sparse leaf
    elseif (Lh.typ == 3)
        Bh = Lh.dat \ Bh;

    % Unknown type    
    else
        error('hmxSolveLower : unavailable case')
    end
  
    
%%% Unavailable     
else
    error('hmxSolveLower.m : unavailable case')
end
end
