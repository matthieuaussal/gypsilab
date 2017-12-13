function Bh = hmxSolveUpper(Bh,Uh)
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
%|    #    |   FILE       : hmxSolveUpper.m                               |
%|    #    |   VERSION    : 0.32                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 25.12.2017                                    |
%| ( === ) |   SYNOPSIS   : Solve upper H-Matrix wit the rule             |
%|  `---'  |                Compr > Full > H-Matrix                       |
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
        
        % Fusion
        Bh = hmxFusion(Bh);
        
    % H-Matrix / Compr -> ---    
    elseif (Bh.typ == 0) && (Uh.typ == 1)   
        error('hmxSolveUpper : unavailable case')
        
    % H-Matrix / Full -> Full   
    elseif (Bh.typ == 0) && (Uh.typ == 2)   
        B      = full(Bh) / Uh.dat;
        Bh     = hmx(Bh.pos{1},Uh.pos{1},Uh.tol);
        Bh.dat = B;
        Bh.typ = 2;
        
        
    % Compr / H-Matrix -> Compr
    elseif (Bh.typ == 1) && (Uh.typ == 0) 
        Bh.dat = {Bh.dat{1} , hmxSolveUpper(Bh.dat{2},Uh)};    
        
    % Compr / Compr -> ---
    elseif (Bh.typ == 1) && (Uh.typ == 1) 
        error('hmxSolveUpper : unavailable case') 
        
    % Compr / Full -> Compr
    elseif (Bh.typ == 1) && (Uh.typ == 2) 
        Bh.dat = {Bh.dat{1} , Bh.dat{2}/Uh.dat}; 
        
        
    % Full / H-Matrix -> Full
    elseif (Bh.typ == 2) && (Uh.typ == 0) 
        Bh.dat = hmxSolveUpper(Bh.dat,Uh);
        
    % Full / Compr -> ---
    elseif (Bh.typ == 2) && (Uh.typ == 1) 
        error('hmxSolveUpper : unavailable case') 
        
    % Full / Full -> Full
    elseif (Bh.typ == 2) && (Uh.typ == 2) 
        Bh.dat = Bh.dat / Uh.dat;    
        
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
    
    % Unknown type    
    else
        error('hmxSolveUpper : unavailable case')
    end
  
    
    
%%% Unavailable
else
    error('hmxSolveUpper.m : unavailable case')
end
end
