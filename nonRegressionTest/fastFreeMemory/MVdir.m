function MV = MVdir(sigma,u,k,beta,tol,V)
%+========================================================================+
%|                                                                        |
%|         OPENFFM - LIBRARY FOR FAST AND FREE MEMORY CONVOLUTION         |
%|           openFfm is part of the GYPSILAB toolbox for Matlab           |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal (c) 2017-2019.                             |
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
%|    #    |   FILE       : MVdir.m                                       |
%|    #    |   VERSION    : 0.6                                           |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2019                                    |
%| ( === ) |   SYNOPSIS   : Non regression test for BEM resolution of     |
%|  `---'  |                a spherical dirichlet scattering problem      |
%+========================================================================+

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [1i*k*beta*S - (Id/2 + D)] = p0 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Matrix-vector product
MV = cell(1,4);

% Loop for each scalar operator ("for" -> "parfor" gives coarse parallelism)
for i = 1:length(MV)
    if (i==1)
        MV{i} = beta/(4*pi) .* integral(sigma,sigma,u,'[exp(ikr)/r]',k,u,tol,V);
    elseif (i==2)
        MV{i} = -1/(4*pi) .* integral(sigma,sigma,u,'grady[exp(ikr)/r]1',k,ntimes(u,1),tol,V);
    elseif (i==3)
        MV{i} = -1/(4*pi) .* integral(sigma,sigma,u,'grady[exp(ikr)/r]2',k,ntimes(u,2),tol,V);
    elseif (i==4)
        MV{i} = -1/(4*pi) .* integral(sigma,sigma,u,'grady[exp(ikr)/r]3',k,ntimes(u,3),tol,V);
    end
end

% Recuperation
MV = cell2mat(MV);
MV = sum(MV,2);
end
