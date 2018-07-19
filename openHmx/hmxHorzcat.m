function Mh = hmxHorzcat(varargin)
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
%|    #    |   FILE       : hmxHorzcat.m                                  |
%|    #    |   VERSION    : 0.42                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 15.07.2018                                    |
%| ( === ) |   SYNOPSIS   : H-Matrix horizontal concatenation             |
%|  `---'  |                                                              |
%+========================================================================+

%%% Input analysis
Ml = varargin{1};
Mr = varargin{2};
if (nargin == 2)
    if isequal(Ml.pos{1},Mr.pos{1})
        Ix = (1:size(Ml,1))';        
        Iy = (1:size(Ml,2)+size(Mr,2))';
    else
        error('hmxHorzcat.m : unavailable case');
    end
else
    Ml = varargin{1};
    Mr = varargin{2};
    Ix = varargin{3};
    Iy = varargin{4};
end

%%% Preparation
% Particles
X = Ml.pos{1}(Ix,:);
Y = [Ml.pos{2} ; Mr.pos{2}];
Y = Y(Iy,:);

% Accuracy
tol = max(Ml.tol,Mr.tol);

% Initialiation
Mh = hmx(X,Y,tol);

% Admissibility
[isfar,Xdim,Ydim] = hmxFar(Mh);

% Y Indices
bool = (Iy<=size(Ml,2));
Iyl  = Iy(bool);
Iyr  = Iy(~bool)-size(Ml,2);

% Compression for far distances
if isfar
    % Compression
    [Al,Bl] = lowrank(Ml,Ix,Iyl);
    [Ar,Br] = lowrank(Mr,Ix,Iyr);
    
    % Concatenation
    A = [Al,Ar];
    B = zeros(size(A,2),length(Iy));
    B(1:size(Al,2),bool)      = Bl;
    B(size(Al,2)+1:end,~bool) = Br;
    
    % Recompression
    [A,B] = hmxQRSVD(A,B,Mh.tol);
    
    % Validation
    flag = (numel(A)+numel(B)) <= prod(size(Mh));
    
else
    flag = 0;
end


%%% Compression
if flag
    % Type
    Mh.typ = 1;
    
    % Low-rank
    Mh.dat = {A,B};


%%%% Full or sparse for smallest box (stopping criterion)
elseif (min(size(Mh)) < 100)
    % Type
    Mh.typ = 2;
    
    % Full
    Mh.dat = zeros(size(Mh));
    
    % Fill data
    Mh.dat(:,bool)  = full(Ml,Ix,Iyl);
    Mh.dat(:,~bool) = full(Mr,Ix,Iyr);


%%% H-Matrix (recursion)
else
    % Type
    Mh.typ = 0;
    
    % Subdivision for X
    [I1,I2] = hmxSubdivide(X,Xdim);
    Mh.row  = {I1 , I1 , I2 , I2};
    
    % Subdivision for Y
    [I1,I2] = hmxSubdivide(Y,Ydim);
    Mh.col  = {I1 , I2 , I1 , I2};
    
    % H-Matrix (recursion)
    for i = 1:4
        % Coordinates
        Ir = Mh.row{i};
        Ic = Mh.col{i};

        % Recursion
        Mh.chd{i} = hmxHorzcat(Ml,Mr,Ix(Ir),Iy(Ic)); 
    end
    
    % Fusion
    Mh = hmxFusion(Mh);    
end
end
