classdef hmx
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
%| or publications in which you use the code. For libFem,     \ \| |/ /   |
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
%| Synopsis   : H-Matrix class definition                                 |
%+========================================================================+
    
properties
    dim = [];         % H-MATRIX DIMENSIONS 
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
        % Empty initialization
        if (length(varargin) == 1)
            Min    = varargin{1};
            Mh.dim = Min.dim;
            Mh.chd = cell(1,4);
            Mh.row = Min.row;
            Mh.col = Min.col;
            Mh.dat = [];
            Mh.typ = [];
            Mh.tol = Min.tol;
        
        % Initialization with dimension and accuracy    
        elseif (length(varargin) == 3)
            Mh.dim = [varargin{1},varargin{2}];
            Mh.chd = cell(1,4);
            Mh.row = cell(1,4);
            Mh.col = cell(1,4);
            Mh.dat = [];
            Mh.typ = [];
            Mh.tol = varargin{3};
       
        % Particles builder with partial abd total pivoting   
        elseif (length(varargin) == 4)
            X     = varargin{1};
            Y     = varargin{2};
            green = varargin{3};
            acc   = varargin{4};
            Mh = hmxBuilder(X,Y,green,acc);
              
        % Compressed builder    
        elseif (length(varargin) == 5)
            X   = varargin{1};
            Y   = varargin{2};
            A   = varargin{3};
            B   = varargin{4};
            acc = varargin{5};
            Mh     = hmx(size(X,1),size(Y,1),acc);
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
            if exist('matlabpool','file') || exist('parpool','file')
                Mh = femHmatrixParallel(Xdof,Ydof,Mx,X,green,Y,My,acc);
            else
                Mh = femHmatrix(Xdof,Ydof,Mx,X,green,Y,My,acc);
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
        if Mh.typ > 0
            B = Mh.dat \ B;
        elseif (Mh.chd{2}.typ == 3)
            B = hmxSolveLower(Mh,B);
        elseif (Mh.chd{3}.typ == 3)
            B = hmxSolveUpper(Mh,B);
        else
            [Lh,Uh] = lu(Mh);
            B = hmxSolveLower(Lh,B);
            B = hmxSolveUpper(Uh,B);
        end
    end
    
    % ITERATIVE SOLVER 
    function B = gmres(Mh,B,restart,tol,maxit)
        B = gmres(@(V) Mh*V,B,restart,tol,maxit);
    end
    
    % STRUCTURE VISUALISATION
    function spy(Mh)
        hmxSpy(Mh);
    end
    
    % DIMENSIONS
    function dim = size(varargin)
        if length(varargin) == 1
            dim = varargin{1}.dim;
        else
            dim = varargin{1}.dim(varargin{2});
        end
    end
    
    % NORM(S)
    
    % MPOWER
end
end
