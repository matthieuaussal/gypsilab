function [Mh,A,B] = hmxSherMorr(Mh)
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
%|    #    |   FILE       : hmxSherMorr.m                                 |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Convert H-Matrix to Shermann Morrison form    |
%|  `---'  |                Sh + A*B                                      |
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
    Mh.dat{1} = zeros(Mh.dim(1),1,class(A)); 
    Mh.dat{2} = zeros(1,Mh.dim(2),class(A));
    
% Full leaf
elseif (Mh.typ == 2)
    A = zeros(Mh.dim(1),0,class(Mh.dat)); 
    B = zeros(0,Mh.dim(2),class(Mh.dat));

% Unknown type
else
    error('hmxSherMorr.m : unavailable case')
end
end
