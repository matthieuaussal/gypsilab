function [Mh,A,B] = hmxSherMorr(Mh)
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
%| Synopsis   : Convert H-Matrix to Shermann Morrison form Sh + A*B       |
%+========================================================================+
    
% H-Matrix (recursion)
if (Mh.typ == 0)
    % Initialization
    A = zeros(Mh.dim(1),0,class(Mh.row{1}));
    B = zeros(0,Mh.dim(2),class(Mh.row{1}));
    l = 0;
    
    % Recursion
    for i = 1:4
        % Children computation
        [Mh.chd{i},Ai,Bi] = hmxSherMorr(Mh.chd{i});
        ri = size(Ai,2);
        
        % Parent incrementation
        A(Mh.row{i},l+1:l+ri) = Ai;
        B(l+1:l+ri,Mh.col{i}) = Bi;
        
        % Rank incrementation
        l = l+ri;
    end
    
    % Recompression
    [A,B] = hmxQRSVD(A,B,Mh.tol);
    
% Compressed leaf
elseif (Mh.typ == 1)
    A = Mh.dat{1}; 
    B = Mh.dat{2};
    Mh.dat{1} = zeros(Mh.dim(1),0,class(A)); 
    Mh.dat{2} = zeros(0,Mh.dim(2),class(A));
    
% Full leaf
elseif (Mh.typ == 2)
    A = zeros(Mh.dim(1),0,class(Mh.dat)); 
    B = zeros(0,Mh.dim(2),class(Mh.dat));

% Sparse leaf
elseif (Mh.typ == 3)
    A = zeros(Mh.dim(1),0); 
    B = zeros(0,Mh.dim(2));

% Unknown type
else
    error('hmxSherMorr.m : unavailable case')
end

end
