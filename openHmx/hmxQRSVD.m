function [A,B] = hmxQRSVD(A,B,tol)
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
%| Synopsis   : QR factorization and SVD compression for low-rank matrices|
%+========================================================================+

% For non-empty matrix
if ~isempty(A)
    % QR Factorisation A = QA * RA;
    [QA,RA] = qr(A,0);
    
    % QR Factorisation Bt = QB * RB
    [QB,RB] = qr(B.',0);
    
    % SVD : U*S*V = RA * RB.'
    [U,S,V] = svd(RA * RB.','econ');
    
    % Indices singular values > tol/10
    ind = find(abs(diag(S)./S(1)) > tol/10);
    
    % Recompression A = QA * U * S
    A = QA * (U(:,ind) * S(ind,ind));
    
    % Recompression B = V' * QB^t
    B = V(:,ind)' * QB.';
end
end
