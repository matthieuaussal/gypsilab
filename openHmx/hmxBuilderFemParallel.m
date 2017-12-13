function [Mh,leaf] = hmxBuilderFemParallel(varargin)
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
%|    #    |   FILE       : hmxBuilderFemParallel.m                       |
%|    #    |   VERSION    : 0.32                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 25.12.2017                                    |
%| ( === ) |   SYNOPSIS   : Parallelisation of finite element builder     |
%|  `---'  |                                                              |
%+========================================================================+

% Fixed input
Xdof  = varargin{1};
Ydof  = varargin{2};
Mx    = varargin{3};
X     = varargin{4};
green = varargin{5};
Y     = varargin{6};
My    = varargin{7};
tol   = varargin{8};

% Variable input
if (length(varargin) == 8)
    leaf = [];
    Nsub = 1;
else
    leaf = varargin{9};
    Nsub = varargin{10};
end

% Initialisation
Mh = hmx(Xdof,Ydof,tol);

% Cores availables for computation
Nlabs = max(1,length(Composite()));

% Pool is 2^n
if (floor(log2(Nlabs)) ~= log2(Nlabs))
    warning('hmxBuilderFemParallel.m : number of pool labs have to be multiple of 2')
    delete(gcp)
    Nlabs = 2^floor(log2(Nlabs));
    parpool(Nlabs);
end

% H-Matrix (recursion)
if (Nlabs >= 2*Nsub)
    % Subdivision for Xdof
    [I1,I2] = hmxSubdivide(Xdof);
    Mh.row  = {I1 , I1 , I2 , I2 };
    
    % Subdivision for Ydof
    [I1,I2] = hmxSubdivide(Ydof);
    Mh.col  = {I1 , I2 , I1 , I2};
    
    % H-Matrix (recursion)
    for i = 1:4
        % Dof indices
        Ir = Mh.row{i};
        Ic = Mh.col{i};
        
        % Fem matrix subdivision and quadratures points indices
        [Mxchd,Ix] = femSubdivideCell(Mx,Ir,'left');
        [Mychd,Iy] = femSubdivideCell(My,Ic,'right');
              
        % Recursion
        [Mh.chd{i},leaf] = hmxBuilderFemParallel(Xdof(Ir,:),Ydof(Ic,:),...
            Mxchd,X(Ix,:),green,Y(Iy,:),Mychd,tol,leaf,2*Nsub);
    end
    
    % Leaves types
    Mh.typ = 0;

% Computation data    
else
    Mh.typ      = -1;
    ind         = size(leaf,1) + 1;
    leaf{ind,1} = Xdof;
    leaf{ind,2} = Ydof;
    leaf{ind,3} = Mx;
    leaf{ind,4} = X;
    leaf{ind,5} = Y;
    leaf{ind,6} = My;
end

% Computation
if Nsub == 1
    % No paralelism
    if Nlabs == 1
         Mh = hmxBuilderFem(leaf{1}, leaf{2}, leaf{3}, leaf{4}, ...
             green,leaf{5},leaf{6},tol);
        
    % Parallelism    
    elseif (size(leaf,1) == Nlabs^2)
        % Distances beetween boxes
        dst = zeros(Nlabs^2,1);
        for i = 1:length(leaf)
        	Xctr   = 0.5 * (min(leaf{i,1},[],1) + max(leaf{i,1},[],1));
            Yctr   = 0.5 * (min(leaf{i,2},[],1) + max(leaf{i,2},[],1));
            dst(i) = sqrt(sum((Yctr-Xctr).^2, 2));
        end

        % Sorted interactions
        [~,I] = sort(dst);
        
        % Loop on workers
        for i = 1:Nlabs
            % Clock
            tps = hmxClock();
            
            % Initialization
            ind = zeros(Nlabs,1);
            tmp = cell(Nlabs,1);
            
            % Parallel computation with optimal
            parfor j = 1:Nlabs
                ind(j) = I((i-1)*Nlabs + j);
                tmp{j} = hmxBuilderFem(leaf{ind(j),1}, leaf{ind(j),2}, ...
                    leaf{ind(j),3}, leaf{ind(j),4}, green, ...
                    leaf{ind(j),5}, leaf{ind(j),6}, tol);
            end
            
            % H-Matrix distribution
            for j = 1:Nlabs
                Mh = hmxDistribution(Mh,tmp{j},ind(j),log2(Nlabs),1);
            end
            
            % Informations
            disp([' ~~> Step ',num2str(i)', ' - Elapsed time is ', ...
                num2str(hmxClock()-tps),' seconds.'])
            
            % Clean composite
            clear tmp;            
        end
    
    % Problem
    else
        error('hmxBuilderFemParallel.m : unavailable case.');
    end
end
end