function Ml = hmxPlus(Ml,Mr)
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
%|    #    |   FILE       : hmxPlus.m                                     |
%|    #    |   VERSION    : 0.31                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 25.11.2017                                    |
%| ( === ) |   SYNOPSIS   : Sum of H-Matrix                               |
%|  `---'  |                                                              |
%+========================================================================+

%%% H-Matrix + H-Matrix -> H-Matrix
if isa(Ml,'hmx') && isa(Mr,'hmx')
    % Check dimensions
    if sum(Ml.dim ~= Mr.dim)
        error('hmxPlus.m : matrix dimensions must agree.')
    end
    
    % H-Matrix + H-Matrix --> H-Matrix (recursion)
    if (Ml.typ == 0) && (Mr.typ == 0)
        % Construction
        for i = 1:4
            Ml.chd{i} = hmxPlus(Ml.chd{i},Mr.chd{i});
        end

        % Fusion
        Ml = hmxFusion(Ml);
        
    else
        % Empty + --- ---> ---
        if issparse(Ml.dat) && isempty(find(Ml.dat,1))
            Ml = Mr;
            
        % --- + Empty ---> ---
        elseif issparse(Mr.dat) && isempty(find(Mr.dat,1))
        
        % H-Matrix + Compr -> H-Matrix
        elseif (Ml.typ == 0) && (Mr.typ == 1)
            Ml = hmxPlusAB(Ml,Mr.dat{1},Mr.dat{2});
             
        % H-Matrix + --- -> ---  
        elseif (Ml.typ == 0) 
            Ml = hmxPlus(Ml,Mr.dat);
    
        % Compr + --- -> --- 
        elseif (Ml.typ == 1)
            Ml = hmxPlusAB(Mr,Ml.dat{1},Ml.dat{2});
            
        % Full + --- -> Full
        elseif (Ml.typ == 2)
            Ml = hmxPlus(Mr,Ml.dat);
            
        % Sparse + --- -> Sparse
        elseif (Ml.typ == 3)
            Ml = hmxPlus(Mr,Ml.dat);
            
        else
            error('hmxPlus : unvailable case')
        end
    end
    
    
%%% H-Matrix + Matrix -> H-Matrix
elseif isa(Ml,'hmx') 
    % Check dimensions
    if sum(Ml.dim ~= size(Mr))
        error('hmxPlus.m : matrix dimensions must agree.')
    end
    
    % H-Matrix (recursion)
    if (Ml.typ == 0)
        % Construction
        for i = 1:4
            Ml.chd{i} = hmxPlus(Ml.chd{i},Mr(Ml.row{i},Ml.col{i}));
        end
        
        % Fusion
        Ml = hmxFusion(Ml);

    % Compressed leaf                             
    elseif (Ml.typ == 1)
        % Save compressed data 
        A = Ml.dat{1};
        B = Ml.dat{2};
        
        % Exchange with matrix
        Ml.dat = Mr;
        if ~issparse(Mr)
            Ml.typ = 2;                         %%%%%%% TYPE CHANGE %%%%%%%
        else
            Ml.typ = 3;                         %%%%%%% TYPE CHANGE %%%%%%%
        end
        
        % Update
        Ml = hmxPlusAB(Ml,A,B);
        
    % Full leaf
    elseif (Ml.typ == 2)
        Ml.dat = Ml.dat + full(Mr);
        
    % Sparse leaf                           
    elseif (Ml.typ == 3)                       
        if issparse(Mr)
            Ml.dat = Ml.dat + Mr;
        else
            Ml.dat = full(Ml.dat) + Mr;
            Ml.typ = 2;                         %%%%%%% TYPE CHANGE %%%%%%%
        end
        
    % Unknown type
    else
        error('hmxPlus.m : unavailable case')
    end

    
%%% Matrix + H-Matrix -> Matrix
elseif isa(Mr,'hmx')
    Ml = hmxPlus(Mr,Ml);
    
%%% Unavailable  
else
    error('hmxPlus.m : unavailable case')
end
end


function Mh = hmxPlusAB(Mh,A,B)
% H-Matrix (recursion)
if (Mh.typ == 0)
    % Construction
    for i = 1:4
        Mh.chd{i} = hmxPlusAB(Mh.chd{i},A(Mh.row{i},:),B(:,Mh.col{i}));
    end
    
    % Fusion
    Mh = hmxFusion(Mh);
    
% Compressed leaf      
elseif (Mh.typ == 1)
    A      = [Mh.dat{1},A];
    B      = [Mh.dat{2};B];
    [A,B]  = hmxQRSVD(A,B,Mh.tol);
    Mh.dat = {A,B};    
    
% Full leaf   
elseif (Mh.typ == 2)
    % Recompression
    [Ah,Bh,flag] = hmxSVD(Mh.dat,Mh.tol);
    if flag && (size(Ah,2) < 0.5*(min(Mh.dim)))
        A      = [Ah,A];
        B      = [Bh;B];
        [A,B]  = hmxQRSVD(A,B,Mh.tol);
        Mh.dat = {A,B};
        Mh.typ = 1;                             %%%%%%% TYPE CHANGE %%%%%%%
    
    % Summation
    else
        Mh.dat = Mh.dat + A*B;
    end
    
% Sparse leaf
elseif (Mh.typ == 3)
    % Empty + Compr -> Compr
    if isempty(Mh.dat)
        Mh.dat = {A,B};
        Mh.typ = 1;                             %%%%%%% TYPE CHANGE %%%%%%%

    % Sparse + Compr -> Compr
    else 
        % Compression
        rk           = min(nnz(Mh.dat)+1,floor(0.5*(min(Mh.dim))));
        [Ah,Bh,flag] = hmxRSVD(Mh.dat,Mh.tol,rk);

        % Update
        if flag
            A      = [A,Ah];
            B      = [B;Bh];
            [A,B]  = hmxQRSVD(A,B,Mh.tol);
            Mh.dat = {A,B};
            Mh.typ = 1;                         %%%%%%% TYPE CHANGE %%%%%%%
        else
            Mh.dat = full(Mh.dat) + A*B;
            Mh.typ = 2;                         %%%%%%% TYPE CHANGE %%%%%%%
        end
    end    
      
% Unknown type
else
    error('hmxPlus.m : unavailable case')
end
end
