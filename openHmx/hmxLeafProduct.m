function MV = hmxLeafProduct(in)
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
%|    #    |   FILE       : hmxLeafProduct.m                              |
%|    #    |   VERSION    : 0.30                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 31.10.2017                                    |
%| ( === ) |   SYNOPSIS   : Parallel matrix-Vector product at leaves      |
%|  `---'  |                                                              |
%+========================================================================+

persistent Mc
persistent N

% Prepare matrix-vector product
if isa(in,'hmx')
    % Leaves blocks
    Ml = hmxLeafOut(in);
    
    % Charge repartition (random)
    Iprm = randperm(length(Ml));
    Ml   = Ml(Iprm);
    
    % Initialization
    Mc = Composite();
    
    % Slicing
    Nloc = ceil(length(Ml)/length(Mc));
    for n = 1:length(Mc)
        ind   = (n-1)*Nloc+1:min(n*Nloc,length(Ml));
        Mc{n} = Ml(ind);
    end
    
    % Dimension
    N = in.dim(1);
    
    % Function handle
    MV = @(V) hmxLeafProduct(V);


% Matrix vector product
else
    % Parallelization
    spmd
        tmp = zeros(N,size(in,2));
        for n = 1:length(Mc)
            I = Mc{n}{1};
            J = Mc{n}{2};
            M = Mc{n}{3};
            if iscell(M)
                tmp(I,:) = tmp(I,:) + M{1} * (M{2} * in(J,:));
            else
                tmp(I,:) = tmp(I,:) + M * in(J,:);
            end
        end
    end
    MV = zeros(N,size(in,2));
    for n = 1:length(tmp)
        MV = MV + tmp{n};
    end
end
end
