function Ml = hmxPlus(Ml,Mr)
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
%| Synopsis   : Matrix summation with H-Matrix                            |
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
    rkMax        = ceil(sqrt(min(Mh.dim))-log(Mh.tol));
    [Ah,Bh,flag] = hmxRSVD(Mh.dat,Mh.tol,rkMax);
     
    % Update
    if flag
        A      = [A,Ah];
        B      = [B;Bh];
        [A,B]  = hmxQRSVD(A,B,Mh.tol);        
        Mh.dat = {A,B};
        Mh.typ = 1;                             %%%%%%% TYPE CHANGE %%%%%%%
        
    else
        Mh.dat = Mh.dat + A*B;
    end
    
% Sparse leaf
elseif (Mh.typ == 3)
    % Indices
    [I,J,V] = find(Mh.dat);
    r       = length(I);
    
    % Empty + Compr -> Compr
    if (r == 0)
        Mh.dat = {A,B};
        Mh.typ = 1;                             %%%%%%% TYPE CHANGE %%%%%%%
        
    % Sparse + Compr -> Compr
    elseif (r > 0) && (r <= 1.1*size(A,2))
        % Low-rank representation
        Ah = zeros(Mh.dim(1),r,class(A));
        Bh = zeros(r,Mh.dim(2),class(B));
        K  = (1:r)';
        
        % Compression
        Ah(sub2ind(size(Ah),I,K)) = ones(size(I),class(A));
        Bh(sub2ind(size(Bh),K,J)) = V;
        
        % Addition
        A = [Ah,A];
        B = [Bh;B];
        
        % Recompression
        [A,B] = hmxQRSVD(A,B,Mh.tol);
        
        % Update
        Mh.dat = {A,B};
        Mh.typ = 1;                             %%%%%%% TYPE CHANGE %%%%%%% 
        
    % Sparse + Compr -> Full
    else
        Mh.dat = full(Mh.dat); 
        Mh.typ = 2;                             %%%%%%% TYPE CHANGE %%%%%%%               
        Mh     = hmxPlusAB(Mh,A,B);
    end

      
% Unknown type
else
    error('hmxPlus.m : unavailable case')
end
end
