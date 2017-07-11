function M = hmxSpy(varargin)
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
%| Synopsis   : Spy H-Matrix architecture                                 |
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
    M       = sparse(Mh.dim(1),Mh.dim(2));
    rk      = size(Mh.dat{1},2);
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
