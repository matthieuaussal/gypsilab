function [Mh,num] = hmxDistribution(Mh,Mp,ind,stp,num)
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
%|    #    |   FILE       : hmxDistribution.m                             |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Distribute H-Matrix parallel leaf to final    |
%|  `---'  |                H-Matrix structure                            |
%+========================================================================+

% H-Matrix (recursion)
if stp && (Mh.typ == 0)
    % Recursion
    for i = 1:4
        [Mh.chd{i},num] = hmxDistribution(Mh.chd{i},Mp,ind,stp-1,num);
    end
    
    % Fusion
    Mh = hmxFusion(Mh);

% Fusioned leaf 
elseif stp
    num = num + 4^stp;
    
% Leaf    
else
    if (num == ind)
        if (Mh.typ == -1)
            Mh = Mp;
        else
            error('hmxDistribution.m : unavailable case.');
        end
    end
    num = num + 1;
end
end