function [Mh,leaf] = hmxBuilderParallel(varargin)
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
%|    #    |   FILE       : hmxBuilderParallel.m                          |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Parallelisation of particle builder           |
%|  `---'  |                                                              |
%+========================================================================+

% Fixed input
X     = varargin{1};
Y     = varargin{2};
green = varargin{3};
tol   = varargin{4};

% Variable input
if (length(varargin) == 4)
    leaf = [];
    Nsub = 1;
else
    leaf = varargin{5};
    Nsub = varargin{6};
end

% Initialisation
Mh = hmx(X,Y,tol);

% Cores availables for computation
Nlabs = max(1,length(Composite()));

% Pool is 2^n
if (floor(log2(Nlabs)) ~= log2(Nlabs))
    warning('hmxBuilderParallel.m : number of pool labs have to be multiple of 2')
    delete(gcp)
    Nlabs = 2^floor(log2(Nlabs));
    parpool(Nlabs);
end

% H-Matrix (recursion)
if (Nlabs >= 2*Nsub)
    % Subdivision for X
    [I1,I2] = hmxSubdivide(X);
    Mh.row  = {I1 , I1 , I2 , I2};
    
    % Subdivision for Y
    [I1,I2] = hmxSubdivide(Y);
    Mh.col  = {I1 , I2 , I1 , I2};
    
    % Single class
    if isa(X,'single')
        for i = 1:4
            Mh.row{i} = single(Mh.row{i});
            Mh.col{i} = single(Mh.col{i});
        end
    end
    
    % H-Matrix (recursion)
    for i = 1:4
        Xi = X(Mh.row{i},:);
        Yi = Y(Mh.col{i},:);
        if isa(green,'function_handle')
            [Mh.chd{i},leaf] = hmxBuilderParallel(Xi,Yi,green,tol,leaf,2*Nsub);
        else
            Mi = green(Mh.row{i},Mh.col{i});
            [Mh.chd{i},leaf] = hmxBuilderParallel(Xi,Yi,Mi,tol,leaf,2*Nsub);
        end
    end
    
    % Leaves types
    Mh.typ = 0;

% Computation data    
else
    Mh.typ      = -1;
    l           = size(leaf,1);
    leaf{l+1,1} = X;
    leaf{l+1,2} = Y;
    leaf{l+1,3} = green;
end

% Computation
if Nsub == 1
    % No parallelism
    if Nlabs == 1
         Mh = hmxBuilder(leaf{1}, leaf{2}, leaf{3}, tol);
        
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
                tmp{j} = hmxBuilder(leaf{ind(j),1}, leaf{ind(j),2},...
                    leaf{ind(j),3}, tol);
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
        error('hmxBuilderParallel.m : unavailable case.');
    end
end
end