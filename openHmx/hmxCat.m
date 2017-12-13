function Mh = hmxCat(dim,Ml,Mr)
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
%|    #    |   FILE       : hmxCat.m                                      |
%|    #    |   VERSION    : 0.32                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 25.12.2017                                    |
%| ( === ) |   SYNOPSIS   : Concatenate H-Matrix along dimensions dim     |
%|  `---'  |                with the rule                                 |
%|         |                H-Matrix > Compr > Full                       |
%+========================================================================+

%%% [H-Matrix & H-Matrix] -> H-Matrix
if ~(isa(Ml,'hmx') && isa(Mr,'hmx'))
    error('hmxConcatenate.m : unavailable case')
elseif (dim == 1) && (size(Ml,2)~=size(Mr,2))
    error('hmxConcatenate.m : Dimensions of matrices being concatenated are not consistent')
elseif (dim == 2) && (size(Ml,1)~=size(Mr,1))
    error('hmxConcatenate.m : Dimensions of matrices being concatenated are not consistent')
elseif (dim > 2)
    error('hmxConcatenate.m : unavailable case')
end
    
% Initialization
if (dim == 1)
    Mh = hmx([Ml.pos{1};Mr.pos{1}],Ml.pos{2},Ml.tol);
elseif (dim == 2)
    Mh = hmx(Ml.pos{1},[Ml.pos{2};Mr.pos{2}],Ml.tol);
else
    error('hmxConcatenate.m : unavailable case')
end


%%% [H-Matrix & H-Matrix] -> H-Matrix
if (Ml.typ == 0) && (Mr.typ == 0)
    % Construction
    for i = 1:4
        Mh.chd{i} = hmxCat(dim,Ml.chd{i},Mr.chd{i});
        if (dim == 1)
            Mh.row{i} = [Ml.row{i};Ml.dim(1)+Mr.row{i}];
            Mh.col{i} = Ml.col{i};
        else
            Mh.row{i} = Ml.row{i};
            Mh.col{i} = [Ml.col{i};Ml.dim(2)+Mr.col{i}];
        end 
    end
    Mh.typ = 0;
    
    % Split for too big full leaf 
    for i = 1:4
        if (Mh.chd{i}.typ == 2) && ~sum(Mh.chd{i}.dim < 200)
            Mh.chd{i} = hmxSplit(dim,Mh.chd{i});
            Mh.chd{i} = hmxFusion(Mh.chd{i});
        end
    end
    
    % Fusion
    Mh = hmxFusion(Mh);
    
    
% [H-Matrix & Compr] -> H-Matrix           
elseif (Ml.typ == 0) && (Mr.typ == 1)
    Mr = hmxSplit(dim,Mr,Ml.row,Ml.col);
    Mh = hmxCat(dim,Ml,Mr);

    
% [H-Matrix & Full] -> H-Matrix 
elseif (Ml.typ == 0) && (Mr.typ == 2)
    Mr = hmxSplit(dim,Mr,Ml.row,Ml.col);
    Mh = hmxCat(dim,Ml,Mr);

    

%%% [Compr & H-matrix] -> H-Matrix
elseif (Ml.typ == 1) && (Mr.typ == 0)
    Ml = hmxSplit(dim,Ml,Mr.row,Mr.col);
    Mh = hmxCat(dim,Ml,Mr);    
  
    
% [Compr & Compr] -> H-matrix or Compr or Full
elseif (Ml.typ == 1) && (Mr.typ == 1)
    % Concatenation
    if (dim == 1)
        A = [Ml.dat{1} , zeros(size(Ml.dat{1},1),size(Mr.dat{1},2)) ;
            zeros(size(Mr.dat{1},1),size(Ml.dat{1},2)) , Mr.dat{1}] ;
        B = [Ml.dat{2} ; Mr.dat{2}];
    else
        A = [Ml.dat{1} , Mr.dat{1}]; 
        B = [Ml.dat{2} , zeros(size(Ml.dat{2},1),size(Mr.dat{2},2)) ;
            zeros(size(Mr.dat{2},1),size(Ml.dat{2},2)) , Mr.dat{2}] ;
    end
    
    % Recompression
    [A,B] = hmxQRSVD(A,B,Ml.tol);
    
    % Update
    Mh.dat = {A,B};
    Mh.typ = 1;
 
    
% [Compr & Full] -> H-matrix or Compr or Full
elseif (Ml.typ == 1) && (Mr.typ == 2)
    if sum(Ml.dim < 200)
        Mh.dat = cat(dim,Ml.dat{1}*Ml.dat{2},Mr.dat);
        Mh.typ = 2;  
    else 
        Mr.dat = {full(Mr.dat),eye(Mr.dim(2))};
        Mr.typ = 1;
        Mh     = hmxCat(dim,Ml,Mr);
    end

    
    
%%% [Full & H-matrix] -> H-Matrix
elseif (Ml.typ == 2) && (Mr.typ == 0)
    Ml = hmxSplit(dim,Ml,Mr.row,Mr.col);
    Mh = hmxCat(dim,Ml,Mr); 
    
        
% [Full & Compr] -> Compr                            
elseif (Ml.typ == 2) && (Mr.typ == 1)       
    if sum(Mr.dim < 200)
        Mh.dat = cat(dim,Ml.dat,Mr.dat{1}*Mr.dat{2});
        Mh.typ = 2;   
    else
        Ml.dat = {full(Ml.dat),eye(Ml.dim(2))};
        Ml.typ = 1;
        Mh     = hmxCat(dim,Ml,Mr);  
        
    end

        
% [Full & Full] -> Full
elseif (Ml.typ == 2) && (Mr.typ == 2)
    Mh.dat = cat(dim,Ml.dat,Mr.dat);
    Mh.typ = 2;   
    
else
    error('hmxCat.m : unavailable case')
end
end



function Ml = hmxSplit(varargin)
% Input
dim = varargin{1};
Ml  = varargin{2};
    
% Subdivision for Xl
if (dim == 1)
    % Row
    [I1,I2] = hmxSubdivide(Ml.pos{1});
    Ml.row  = {I1 , I1 , I2 , I2};
    
    % Column
    if (nargin == 4)
        Ml.col = varargin{4};
    else
        [I1,I2] = hmxSubdivide(Ml.pos{2});
        Ml.col  = {I1 , I2 , I1 , I2};
    end
    
% Subdivision for Yl
else
    % Row
    if (nargin == 4)
        Ml.row = varargin{3};
    else
        [I1,I2] = hmxSubdivide(Ml.pos{1});
        Ml.row  = {I1 , I1 , I2 , I2};
    end
    
    % Column
    [I1,I2] = hmxSubdivide(Ml.pos{2});
    Ml.col  = {I1 , I2 , I1 , I2};
end

% Create H-Matrix
for i = 1:4
    % Initialize
    Ml.chd{i} = hmx(Ml.pos{1}(Ml.row{i},:),Ml.pos{2}(Ml.col{i},:),Ml.tol);
    
    % H-Matrix
    if (Ml.typ == 0)
        error('hmxSplit.m : unavailable case')

    % Compressed leaf
    elseif (Ml.typ == 1)
        % Subdivision
        A = Ml.dat{1}(Ml.row{i},:);
        B = Ml.dat{2}(:,Ml.col{i});
        
        % Recompression
        [A,B] = hmxQRSVD(A,B,Ml.tol);
        
        % Update
        Ml.chd{i}.dat = {A,B};
        Ml.chd{i}.typ = 1;
        
    % Full leaf    
    elseif (Ml.typ == 2)        
        Ml.chd{i}.dat = Ml.dat(Ml.row{i},Ml.col{i});
        Ml.chd{i}.typ = 2;   
        
    else
        error('hmxSplit.m : unavailable case')
    end
end

% Nullify leaf
Ml.dat = [];
Ml.typ = 0;
end
