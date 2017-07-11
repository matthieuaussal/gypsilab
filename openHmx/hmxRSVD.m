function [A,B,flag] = hmxRSVD(M,tol,rk)
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
%|              Antoine Liutkus - Inria                                   |
%| Creation   : 14.03.17                                                  |
%| Last modif : 21.06.17                                                  |
%| Synopsis   : RSVD from Halko et al. 'finding structure with randomness |
%+========================================================================+

% Dimensions
[m,n] = size(M);
p     = min(2*rk,n);

% Randomized
X  = randn(n,p,class(M));
Y  = M*X;
W1 = orth(Y);
B  = W1'*M;

% Truncated SVD
[W2,S,V] = svd(B,'econ');
U        = W1*W2;

% SVD final matrix
rk = min(rk,size(U,2));
U  = U(:,1:rk);
S  = S(1:rk,1:rk);
V  = V(:,1:rk);

% No values
if numel(S) == 0
   A = zeros(m,0);
   B = zeros(0,n);
   flag = 1;
   return
end

% Final low-rank with fixed accuracy
I = find(abs(diag(S))/S(1) >= tol/10);
if (length(I) == rk)
    A    = [];
    B    = [];
    flag = 0;
else
    A    = U(:,I);
    B    = S(I,I) * V(:,I)';
    flag = 1;
end

end
