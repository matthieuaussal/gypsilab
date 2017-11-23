function [Mh,leaf] = femHmxBuilderParallel(varargin)
%+========================================================================+
%|                                                                        |
%|              OPENFEM - LIBRARY FOR FINITE ELEMENT METHOD               |
%|           openFem is part of the GYPSYLAB toolbox for Matlab           |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal & Francois Alouges (c) 2015-2017           |
%| PROPERTY  : Centre de Mathematiques Appliquees, Ecole polytechnique,   |
%| route de Saclay, 91128 Palaiseau, France. All rights reserved          |
%| LICENCE   :                                                            |
%| This program is a full open source software. It is distributed in the  |
%| hope that it will be useful, but without any warranty. This program is |
%| free only for public purpose: you can modify and/or adapt it only if   |
%| you share all your contributions and results, in open source and under |
%| the same licence. For other uses, you may contact us to obtain a       |
%| suitable license. Please acknowledge the gypsilab toolbox in programs  |
%| or publications in which you use it.                                   |
%| CONTACT   : matthieu.aussal@polytechnique.edu                          |
%|             francois.alouges@polytechnique.edu                         |
%| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
%|________________________________________________________________________|
%|   '&`   |                                                              |
%|    #    |   FILE       : femHmxBuilderParallel.m                       |
%|    #    |   VERSION    : 0.31                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 25.11.2017                                    |
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
Mh = hmx(size(Xdof,1),size(Ydof,1),tol);

% Cores availables for computation
Nlabs = max(1,length(Composite()));

% Pool is 2^n
if (floor(log2(Nlabs)) ~= log2(Nlabs))
    warning('femHmxBuilderParallel.m : number of pool labs have to be multiple of 2')
    delete(gcp)
    Nlabs = 2^floor(log2(Nlabs));
    parpool(Nlabs);
end

% H-Matrix (recursion)
if (Nlabs >= 2*Nsub)
     % Subdivision for Xdof
    ind    = hmxSubdivide(Xdof);
    Mh.row = {find(ind),find(ind),find(~ind),find(~ind)};
    
    % Subdivision for Ydof
    ind    = hmxSubdivide(Ydof);
    Mh.col = {find(ind),find(~ind),find(ind),find(~ind)};
    
    % Single class
    if isa(Xdof,'single')
        for i = 1:4
            Mh.row{i} = single(Mh.row{i});
            Mh.col{i} = single(Mh.col{i});
        end
    end
    
    % H-Matrix (recursion)
    for i = 1:4
        % Dof indices
        Ir = Mh.row{i};
        Ic = Mh.col{i};
        
        % Fem matrix subdivision and quadratures points indices
        [Mxchd,Ix] = femSubdivideCell(Mx,Ir,'left');
        [Mychd,Iy] = femSubdivideCell(My,Ic,'right');
              
        % Recursion
        [Mh.chd{i},leaf] = femHmxBuilderParallel(Xdof(Ir,:),Ydof(Ic,:),...
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
         Mh = femHmatrix(leaf{1}, leaf{2}, leaf{3}, leaf{4}, ...
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
                tmp{j} = femHmxBuilder(leaf{ind(j),1}, leaf{ind(j),2}, ...
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
        error('femHmxBuilderParallel.m : unavailable case.');
    end
end
end