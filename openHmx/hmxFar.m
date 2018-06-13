function [bool,Xdim,Ydim] = hmxFar(Mh)
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
%|    #    |   FILE       : hmxFar.m                                      |
%|    #    |   VERSION    : 0.41                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 01.04.2018                                    |
%| ( === ) |   SYNOPSIS   : H-Matrix far boolean                          |
%|  `---'  |                                                              |
%+========================================================================+

% Particles box X
X        = Mh.pos{1};
Xmin     = min(X,[],1);
Xmax     = max(X,[],1);
Xctr     = 0.5*(Xmin+Xmax);
Xdgl     = Xmax-Xmin;
[~,Xdim] = max(Xdgl);

% Particles box Y
Y        = Mh.pos{2};
Ymin     = min(Y,[],1);
Ymax     = max(Y,[],1);
Yctr     = 0.5*(Ymin+Ymax);
Ydgl     = Ymax-Ymin;
[~,Ydim] = max(Ydgl);

% Separated particles set
bool = sum( (abs(Yctr-Xctr)>=0.75*(Xdgl+Ydgl)) & (Xdgl>0) & (Ydgl > 0));
end
