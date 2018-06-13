function I = domIntegral8(data)
%+========================================================================+
%|                                                                        |
%|              OPENDOM - LIBRARY FOR NUMERICAL INTEGRATION               |
%|           openDom is part of the GYPSILAB toolbox for Matlab           |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal & Francois Alouges (c) 2017-2018.          |
%| PROPERTY  : Centre de Mathematiques Appliquees, Ecole polytechnique,   |
%| route de Saclay, 91128 Palaiseau, France. All rights reserved.         |
%| LICENCE   : This program is free software, distributed in the hope that|
%| it will be useful, but WITHOUT ANY WARRANTY. Natively, you can use,    |
%| redistribute and/or modify it under the terms of the GNU General Public|
%| License, as published by the Free Software Foundation (version 3 or    |
%| later,  http://www.gnu.org/licenses). For private use, dual licencing  |
%| is available, please contact us to activate a "pay for remove" option. |
%| CONTACT   : matthieu.aussal@polytechnique.edu                          |
%|             francois.alouges@polytechnique.edu                         |
%| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
%|                                                                        |
%| Please acknowledge the gypsilab toolbox in programs or publications in |
%| which you use it.                                                      |
%|________________________________________________________________________|
%|   '&`   |                                                              |
%|    #    |   FILE       : domIntegral8.m                                |
%|    #    |   VERSION    : 0.41                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 01.04.2018                                    |
%| ( === ) |   SYNOPSIS   : Numerical integation with 8 input arguments   |
%|  `---'  |                                                              |
%+========================================================================+

%%% FAST AND FURIOUS METHOD WITH BOUNDARY ELEMENT OPERATOR  --> \int_{mesh(x)} \int_{mesh(y)} psi(x)' f(x,y) psi(y) dxdy
% Domain with quadrature
Xdom   = data{1};
[X,Wx] = Xdom.qud;
Nx     = size(X,1);
Wx     = spdiags(Wx,0,Nx,Nx);

% Domain with quadrature
Ydom   = data{2};
[Y,Wy] = Ydom.qud;
Ny     = size(Y,1);
Wy     = spdiags(Wy,0,Ny,Ny);

% Finite element matrix with integration
u  = data{3};
Mu = u.uqm(Xdom);
if iscell(Mu)
    Mu{1} = Mu{1}' * Wx;
    Mu{2} = Mu{2}' * Wx;
    Mu{3} = Mu{3}' * Wx;
else
    Mu = Mu' * Wx;
end

% Green kernel
green = data{4};

% Wave number
k = data{5};

% Finite element matrix with integration
v  = data{6};
Mv = v.uqm(Ydom);
if iscell(Mv)
    Mv{1} = Wy * Mv{1};
    Mv{2} = Wy * Mv{2};
    Mv{3} = Wy * Mv{3};
else
    Mv = Wy * Mv;
end

% Accuracy
tol = data{7};

% Block size
V = data{8};

% FFM matrix-vector product
if ~iscell(Mu) && ~iscell(green) && ~iscell(Mv)
    I = Mu * ffmProduct(X,Y,Mv*V,green,k,tol);
    
elseif iscell(Mu) && ~iscell(green) && iscell(Mv)
    I = 0;
    for i = 1:3
        I = I + Mu{i} * ffmProduct(X,Y,Mv{i}*V,green,k,tol);
    end
    
elseif iscell(Mu) && iscell(green) && ~iscell(Mv)
    I = 0;
    for i = 1:3
        I = I + Mu{i} * ffmProduct(X,Y,Mv*V,green{i},k,tol);
    end
    
elseif ~iscell(Mu) && iscell(green) && iscell(Mv)
    I = 0;
    for i = 1:3
        I = I + Mu * ffmProduct(X,Y,Mv{i}*V,green{i},k,tol);
    end
    
elseif iscell(Mu) && iscell(green) && iscell(Mv)
    I   = 0;
    ind = [1 2 3 ; 1 3 2 ; 2 3 1 ; 2 1 3 ; 3 1 2 ; 3 2 1];
    sgn = [+1 -1 +1 -1 +1 -1];
    for i = 1:6
        I = I + sgn(i) .* (Mu{ind(i,1)} * ...
            ffmProduct(X,Y,Mv{ind(i,3)}*V,green{ind(i,2)},k,tol));
    end
    
else
    error('domIntegral8.m : unavailable case')
end
end
