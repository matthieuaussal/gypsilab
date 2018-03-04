classdef hmx
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
%|    #    |   FILE       : hmx.m                                         |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : H-Matrix class definition and functions       |
%|  `---'  |                                                              |
%+========================================================================+

properties
    dim = [];         % H-MATRIX DIMENSIONS 
    pos = [];         % COORDINATES POSITIONS (X,Y)
    chd = [];         % CHILDREN (M11 M12 M21 M22)
    row = [];         % CHILDREN ROWS INDICES
    col = [];         % CHILDREN COLUMNS INDICES
    dat = [];         % LEAF DATA
    typ = [];         % LEAF TYPE (0=H-MATRIX ; 1=COMPRESSED ; 2=FULL) 
    tol = [];         % COMPRESSORS ACCURACY
end

methods
    % CONSTRUCTOR
    function Mh = hmx(varargin)
        % Initialization with dimension and accuracy    
        if (length(varargin) == 3)
            Mh.dim = [length(varargin{1}),length(varargin{2})];
            Mh.pos = {varargin{1},varargin{2}};
            Mh.chd = cell(1,4);
            Mh.row = cell(1,4);
            Mh.col = cell(1,4);
            Mh.dat = [];
            Mh.typ = [];
            Mh.tol = varargin{3};
       
        % Particles builder with partial and total pivoting   
        elseif (length(varargin) == 4)
            X     = varargin{1};
            Y     = varargin{2};
            green = varargin{3};
            acc   = varargin{4};
            if exist('parpool','file')
                Mh = hmxBuilderParallel(X,Y,green,acc);
            else
                Mh = hmxBuilder(X,Y,green,acc);
            end
            
        % Compressed builder    
        elseif (length(varargin) == 5)
            X      = varargin{1};
            Y      = varargin{2};
            A      = varargin{3};
            B      = varargin{4};
            acc    = varargin{5};
            Mh     = hmx(X,Y,acc);
            Mh.dat = {A,B};
            Mh.typ = 1;
            
       % Finite element builder    
        elseif (length(varargin) == 8)
            Xdof  = varargin{1};
            Ydof  = varargin{2};
            Mx    = varargin{3};
            X     = varargin{4};
            green = varargin{5};
            Y     = varargin{6};
            My    = varargin{7};
            acc   = varargin{8};
            if exist('parpool','file')
                Mh = hmxBuilderFemParallel(Xdof,Ydof,Mx,X,green,Y,My,acc);
            else
                Mh = hmxBuilderFem(Xdof,Ydof,Mx,X,green,Y,My,acc);
            end
            
        else
            error('hmx.m : undefined constructor case')
        end
    end    
    
    % FULL CONVERSION
    function M = full(Mh)
        M = hmxFull(Mh);
    end
    
    % SPARSE CONVERSION
    function M = sparse(Mh)
        M = hmxSparse(Mh);
    end

    % SINGLE CONVERSION
    function M = single(Mh)
        M = hmxSingle(Mh);
    end

    % DOUBLE CONVERSION
    function M = double(Mh)
        M = hmxDouble(Mh);
    end

    % TRANSPOSITION
    function Mh = transpose(Mh)
        Mh = hmxTranspose(Mh);
    end
    
    % TRANSPOSITION CONJUGATE
    function Mh = ctranspose(Mh)
        Mh = hmxCtranspose(Mh);
    end
    
    % SCALAR PRODUCT
    function Mh = times(Ml,Mr)
        Mh = hmxTimes(Ml,Mr);
    end
    
    % UMINUS 
    function Mh = uminus(Mh)
        Mh = (-1).*Mh;
    end
    
    % ADDITION
    function Mh = plus(Ml,Mr)
        Mh = hmxPlus(Ml,Mr);
    end
    
    % SUBSTRACTION
    function Mh = minus(Ml,Mr)
       Mh = Ml + (-Mr); 
    end
    
    % MATRIX PRODUCT
    function Mh = mtimes(Ml,Mr)
        if isnumeric(Ml) && (numel(Ml) == 1)
            Mh = Ml .* Mr;
        elseif isnumeric(Mr) && (numel(Mr) == 1)
            Mh = Ml .* Mr;
        else
            Mh = hmxMtimes(Ml,Mr);
        end
    end
    
    % INVERSION
    function Mh = inv(Mh)
        Mh = hmxInv(Mh);
    end
    
    % CHOLESKY FACTORISATION
    function Mh = chol(Mh)
        Mh = hmxChol(Mh);
    end
    
    % LDLt FACTORISATION
    function [Mh,Dh] = ldl(Mh)
        [Mh,Dh] = hmxLdl(Mh);
    end
    
    % LU FACTORISATION
    function [Lh,Uh] = lu(Mh)
        [Lh,Uh] = hmxLU(Mh);
    end
    
    % EXACT SOLVER 
    function B = mldivide(Mh,B)
        if (Mh.typ == 2)
            B = Mh.dat \ B;
        elseif (Mh.chd{2}.typ == 1) && isempty(Mh.chd{2}.dat{1})
            B = hmxSolveLower(Mh,B);
        elseif (Mh.chd{3}.typ == 1) && isempty(Mh.chd{3}.dat{1})
            B = hmxSolveUpper(Mh,B);
        else
            [Lh,Uh] = lu(Mh);
            B = hmxSolveLower(Lh,B);
            B = hmxSolveUpper(Uh,B);
        end
    end
    
    % STRUCTURE VISUALISATION
    function spy(Mh)
        hmxSpy(Mh);
    end
    
    % PLOT POSITIONS
    function plot3(Mh)
        plot3(Mh.pos{1}(:,1),Mh.pos{1}(:,2),Mh.pos{1}(:,3),'bo',...
            Mh.pos{2}(:,1),Mh.pos{2}(:,2),Mh.pos{2}(:,3),'*r');
    end
    
    % DIMENSIONS
    function dim = size(varargin)
        if length(varargin) == 1
            dim = varargin{1}.dim;
        else
            dim = varargin{1}.dim(varargin{2});
        end
    end
    
    % VERTICAL CONCATENATION
    function Mh = vertcat(varargin)
        Mh = varargin{1};
        for i = 2:nargin
            Mh = cat(1,Mh,varargin{i});
        end
    end

    % HORIZONTAL CONCATENATION
    function Mh = horzcat(varargin)
        Mh = varargin{1};
        for i = 2:nargin
            Mh = cat(2,Mh,varargin{i});
        end
    end
    
    % CONCATENATION
    function Mh = cat(dim,Ml,Mr)
        Mh = hmxCat(dim,Ml,Mr);
    end
    
    % NORM(S)
    
    % MPOWER
end
end
