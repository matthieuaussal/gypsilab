function M = hmxSpy(varargin)
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
%|    #    |   FILE       : hmxSpy.m                                      |
%|    #    |   VERSION    : 0.30                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 31.10.2017                                    |
%| ( === ) |   SYNOPSIS   : Spy H-Matrix architecture                     |
%|  `---'  |                                                              |
%+========================================================================+

% Input data
Mh = varargin{1};

% H-Matrix (recursion)
if (Mh.typ == 0)
    A = hmxSpy(Mh.chd{1},[]);
    B = hmxSpy(Mh.chd{2},[]);
    C = hmxSpy(Mh.chd{3},[]);
    D = hmxSpy(Mh.chd{4},[]);
    M = [A,B;C,D];

% Compressed leaf
elseif (Mh.typ == 1)
    M = sparse(Mh.dim(1),Mh.dim(2));
    if iscell(Mh.dat)
        rk = max(1,size(Mh.dat{1},2));
    else
        rk = 1;
    end
    M(:,1)  = 1;
    M(:,rk) = 1;
    M(1,:)  = 1;
    M(rk,:) = 1;
       
% Full leaf
elseif (Mh.typ == 2)
    M        = sparse(Mh.dim(1),Mh.dim(2));
    M(:,1)   = 2;
    M(:,end) = 2;
    M(1,:)   = 2;
    M(end,:) = 2;
    
% Sparse leaf
elseif (Mh.typ == 3)
    M = 3*(abs(Mh.dat)>0);
    if isempty(find(M,1))
        M(:,1)   = 4;
        M(:,end) = 4;
        M(1,:)   = 4;
        M(end,:) = 4;
    end

% Unknown type    
else
    error('hmxSpy.m : unavailable case')
end

% Graphical representation
if (length(varargin) == 1)
    spy(M==1,'b')
    hold on
    spy(M==2,'r')
    spy(M==3,'m')
    spy(M==4,'g')
    hold off
    title('openHmx : H-Matrix structure');
end
end
