function Bh = hmxSolveLower(Lh,Bh)
%+========================================================================+
%|                                                                        |
%|         OPENHMX - LIBRARY FOR H-MATRIX COMPRESSION AND ALGEBRA         |
%|           openHmx is part of the GYPSILAB toolbox for Matlab           |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal (c) 2017-2018.                             |
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
%|    #    |   FILE       : hmxSolveLower.m                               |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Solve lower H-Matrix with the rule            |
%|  `---'  |                Compr > Full > H-Matrix                       |
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
        
        % Fusion
        Bh = hmxFusion(Bh);
        
    % H-Matrix \ Compr -> Compr  
    elseif (Lh.typ == 0) && (Bh.typ == 1) 
        Bh.dat = { hmxSolveLower(Lh,Bh.dat{1}) , Bh.dat{2} };
        
    % H-Matrix \ Full -> Full
    elseif (Lh.typ == 0) && (Bh.typ == 2) 
        Bh.dat = hmxSolveLower(Lh,Bh.dat);
        
     % Compr \ --- -> ---
    elseif (Lh.typ == 1) 
        error('hmxSolveLower : unavailable case')
                
        
    % Full \ H-Matrix -> Full
    elseif (Lh.typ == 2) && (Bh.typ == 0)
        B      = Lh.dat \ full(Bh);
        Bh     = hmx(Lh.pos{2},Bh.pos{2},Lh.tol);
        Bh.dat = B;
        Bh.typ = 2;
        
    % Full \ Compr -> Compr
    elseif (Lh.typ == 2) && (Bh.typ == 1)
        Bh.dat = { Lh.dat\Bh.dat{1} , Bh.dat{2} };
        
    % Full \ Full -> Full
    elseif (Bh.typ == 2)
        Bh.dat = Lh.dat \ Bh.dat;
        
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

    % Unknown type    
    else
        error('hmxSolveLower : unavailable case')
    end
  
    
%%% Unavailable     
else
    error('hmxSolveLower.m : unavailable case')
end
end
