function MV = MVneu(sigma,u,k,beta,tol,V)
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
%|    #    |   FILE       : MVneu.m                                       |
%|    #    |   VERSION    : 0.6                                           |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2019                                    |
%| ( === ) |   SYNOPSIS   : Non regression test for BEM resolution of     |
%|  `---'  |                a spherical neumann scattering problem        |
%+========================================================================+

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [1i*k*beta*(-Id/2 + Dt) - H] = dnp0 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Matrix-vector product
MV = cell(1,9);

% Loop for each scalar operator ("for" -> "parfor" gives coarse parallelism)
for i = 1:length(MV)
    if (i==1)
        MV{i} = -k^2/(4*pi) .* integral(sigma,sigma,ntimes(u,1),'[exp(ikr)/r]',k,ntimes(u,1),tol,V);
    elseif (i==2)
        MV{i} = -k^2/(4*pi) .* integral(sigma,sigma,ntimes(u,2),'[exp(ikr)/r]',k,ntimes(u,2),tol,V);
    elseif (i==3)
        MV{i} = -k^2/(4*pi) .* integral(sigma,sigma,ntimes(u,3),'[exp(ikr)/r]',k,ntimes(u,3),tol,V);
    elseif (i==4)
        MV{i} = 1/(4*pi) .* integral(sigma,sigma,nxgrad(u,1),'[exp(ikr)/r]',k,nxgrad(u,1),tol,V);
    elseif (i==5)
        MV{i} = 1/(4*pi) .* integral(sigma,sigma,nxgrad(u,2),'[exp(ikr)/r]',k,nxgrad(u,2),tol,V);
    elseif (i==6)
        MV{i} = 1/(4*pi) .* integral(sigma,sigma,nxgrad(u,3),'[exp(ikr)/r]',k,nxgrad(u,3),tol,V);
    elseif (i==7)
        MV{i} = beta/(4*pi) .* integral(sigma,sigma,ntimes(u,1),'gradx[exp(ikr)/r]1',k,u,tol,V);
    elseif (i==8)
        MV{i} = beta/(4*pi) .* integral(sigma,sigma,ntimes(u,2),'gradx[exp(ikr)/r]2',k,u,tol,V);
    elseif (i==9)
        MV{i} = beta/(4*pi) .* integral(sigma,sigma,ntimes(u,3),'gradx[exp(ikr)/r]3',k,u,tol,V);
    end
end

% Recuperation
MV = cell2mat(MV);
MV = sum(MV,2);
end
