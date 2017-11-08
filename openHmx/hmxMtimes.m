function Mh = hmxMtimes(Ml,Mr)
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
%|    #    |   FILE       : hmxMtimes.m                                   |
%|    #    |   VERSION    : 0.30                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 31.10.2017                                    |
%| ( === ) |   SYNOPSIS   : Product of H-Matrix                           |
%|  `---'  |                                                              |
%+========================================================================+

%%% H-Matrix * H-Matrix --> H-Matrix
if isa(Ml,'hmx') && isa(Mr,'hmx')
    % Check dimensions
    if (Ml.dim(2) ~= Mr.dim(1))
        error('hmxMtimes.m : matrix dimensions must agree.')
    end
        
    % Initialisation
    Mh = hmx(Ml.dim(1),Mr.dim(2),Ml.tol);
    
    % (H-Matrix * H-Matrix) --> H-Matrix   (recursion)
    if (Ml.typ==0) && (Mr.typ==0)
        % Bloc product indices
        I = [1 2; 1 2; 3 4; 3 4];
        J = [1 3; 2 4; 1 3; 2 4];
        
        % Bloc product
        for i = 1:4
            Mh.chd{i} = hmxMtimes(Ml.chd{I(i,1)},Mr.chd{J(i,1)}) + ...
                hmxMtimes(Ml.chd{I(i,2)},Mr.chd{J(i,2)});
            Mh.row{i} = Ml.row{I(i,1)};
            Mh.col{i} = Mr.col{J(i,2)};
        end
        Mh.typ = 0;
        
        % Fusion
        Mh = hmxFusion(Mh);
        
    else
        % Empty * --- ---> ---
        if issparse(Ml.dat) && isempty(find(Ml.dat,1))
            Mh.dat = sparse(Mh.dim(1),Mh.dim(2));
            Mh.typ = 3;
            
        % --- * Empty ---> ---
        elseif issparse(Mr.dat) && isempty(find(Mr.dat,1))
            Mh.dat = sparse(Mh.dim(1),Mh.dim(2));
            Mh.typ = 3;
        
        % H-Matrix * Compr --> Compr
        elseif (Ml.typ==0) && (Mr.typ==1)
            Mh.dat = {hmxMtimes(Ml,Mr.dat{1}) , Mr.dat{2}};
            Mh.typ = 1;
            
        % H-Matrix * --- --> ---
        elseif (Ml.typ==0)
            Mh.dat = hmxMtimes(Ml,Mr.dat);
            if ~issparse(Mh.dat)
                Mh.typ = 2;
            else
                Mh.typ = 3;
            end    
            
        % Compr * --- --> ---
        elseif (Ml.typ==1)            
            % Compr * H-Matrix --> Compr
            if (Mr.typ==0)
                Mh.dat = {Ml.dat{1} , hmxMtimes(Ml.dat{2},Mr)};
                Mh.typ = 1;
                
            % Compr * Compr --> Compr
            elseif (Mr.typ==1)
                Mh.dat = {Ml.dat{1} , (Ml.dat{2} * Mr.dat{1}) * Mr.dat{2}};
                Mh.typ = 1;
                
            % Compr * Full --> Compr
            elseif (Mr.typ==2)
                Mh.dat = {Ml.dat{1} , Ml.dat{2} * Mr.dat};
                Mh.typ = 1;
                
            % Compr * Sparse --> Compr
            elseif (Mr.typ==3)
                Mh.dat = {Ml.dat{1} , Ml.dat{2} * Mr.dat};
                Mh.typ = 1;
                
            else
                error('hmxMtimes : unvailable case')
            end
            
        % Full * Compr --> compr
        elseif (Ml.typ==2) && (Mr.typ==1)
            Mh.dat = {Ml.dat*Mr.dat{1} , Mr.dat{2}};
            Mh.typ = 1;
            
        % Full * --- --> ---
        elseif (Ml.typ==2)
            Mh.dat = hmxMtimes(Ml.dat,Mr);
            Mh.typ = 2;    
            
        % Sparse * Compr --> compr
        elseif (Ml.typ==3) && (Mr.typ==1)
            Mh.dat = {Ml.dat*Mr.dat{1} , Mr.dat{2}};
            Mh.typ = 1;    
            
        % Sparse * --- --> ---
        elseif (Ml.typ==3)
            Mh.dat = hmxMtimes(Ml.dat,Mr);
            if ~issparse(Mh.dat)
                Mh.typ = 2;
            else
                Mh.typ = 3;
            end 
            
        else
            error('hmxMtimes : unvailable case')
        end
    end
    
    
%%% H-Matrix * Matrix --> Matrix
elseif isa(Ml,'hmx')
    % Check dimensions
    if (Ml.dim(2) ~= size(Mr,1))
        error('hmxMtimes.m : matrix dimensions must agree.')
    end
    
    % H-Matrix (recursion)
    if (Ml.typ == 0)
        % Recursion
        tmp = cell(1,4);
        spr = 0;
        for i = 1:4
            tmp{i} = hmxMtimes(Ml.chd{i},Mr(Ml.col{i},:));
            spr    = spr + issparse(tmp{i});
        end
        
        % Initialization
        if (spr==4)
            Mh = sparse(Ml.dim(1),size(Mr,2));
        else
            Mh = zeros(Ml.dim(1),size(Mr,2),class(Mr));
        end
        
        % Update
        for i = 1:4
            Mh(Ml.row{i},:) = Mh(Ml.row{i},:) + tmp{i};
        end
        
   % Compressed leaf
    elseif (Ml.typ == 1)
        Mh = Ml.dat{1} * (Ml.dat{2} * Mr);
        
    % Full leaf
    elseif (Ml.typ == 2)
        Mh = Ml.dat * full(Mr);

    % Sparse leaf
    elseif (Ml.typ == 3)
        if issparse(Mr)
            Mh = Ml.dat * Mr;
        else
            Mh = full(Ml.dat) * Mr;
         end
        
    % Unknown type
    else
        error('hmxMtimes.m : unavailable case')
    end
    
    
%%% Matrix * H-Matrix --> Matrix
elseif isa(Mr,'hmx')
    Mh = hmxMtimes(Mr.',Ml.').';
    
    
%%% Unavailable    
else
    error('hmxMtimes.m : unavailable case')
end
end
