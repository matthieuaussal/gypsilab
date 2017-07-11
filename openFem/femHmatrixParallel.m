function [Mh,leaf] = femHmatrixParallel(varargin)
%+========================================================================+
%|                                                                        |
%|                  OPENFEM, FINITE AND BOUNDARY ELEMENT                  |
%|              openFem is part of GYPSYLAB toolbox - v0.20               |
%|                                                                        |
%| Copyright (c) 20015-2017, Ecole polytechnique, all rights reserved.    |
%| Licence Creative Commons BY-NC-SA 4.0, Attribution, NonCommercial and  |
%| ShareAlike (see http://creativecommons.org/licenses/by-nc-sa/4.0/).    |
%| This software is the property from Centre de Mathematiques Appliquees  |
%| de l'Ecole polytechnique, route de Saclay, 91128 Palaiseau, France.    |
%|                                                            _   _   _   |
%| Please acknowledge the GYPSILAB toolbox in programs       | | | | | |  |
%| or publications in which you use the code. For openFem,    \ \| |/ /   |
%| we suggest as reference :                                   \ | | /    |
%| [1] : www.cmap.polytechnique.fr/~aussal/gypsilab             \   /     |
%|                                                               | |      |
%|_______________________________________________________________|_|______|
%| Author(s)  : Matthieu Aussal - CMAP, Ecole polytechnique               |
%| Creation   : 14.03.17                                                  |
%| Last modif : 21.06.17                                                  |
%| Synopsis   : Paralelization of H-matrix finite element builder         |
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
Mh    = hmx(size(Xdof,1),size(Ydof,1),tol);
Ileaf = Composite();

% Cores availables for computation
Nlabs = max(1,length(Ileaf));

% Security check
if (floor(log2(Nlabs)) ~= log2(Nlabs))
    error('femHmatrixParallel.m : number of pool labs have to be multiple of 2')
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
        
        % Quadrature indices
        if iscell(Mx) && iscell(My)
            Vx = 0;
            Vy = 0;
            for j = 1:3
                Vx = Vx + Mx{j}(:,Ir) * (1+rand(length(Ir),1));
                Vy = Vy + My{j}(:,Ic) * (1+rand(length(Ic),1));
            end
        else
           Vx = Mx(:,Ir) * (1+rand(length(Ir),1));
           Vy = My(:,Ic) * (1+rand(length(Ic),1));
        end
        Ix = find(Vx);
        Iy = find(Vy);
        
        % Matrix subdivision
        if iscell(Mx) && iscell(My)
            Mxchd = cell(1,3);
            Mychd = cell(1,3);
            for j = 1:3
                Mxchd{j} = Mx{j}(Ix,Ir);
                Mychd{j} = My{j}(Iy,Ic);
            end
        else
            Mxchd = Mx(Ix,Ir);
            Mychd = My(Iy,Ic);
        end

        % Recursion
        [Mh.chd{i},leaf] = femHmatrixParallel(Xdof(Ir,:),Ydof(Ic,:),...
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
            
            % Optimal repartition
            spmd
                ind = I((i-1)*Nlabs + labindex);
                tmp = femHmatrix(leaf{ind,1}, leaf{ind,2}, ...
                    leaf{ind,3}, leaf{ind,4}, green, ...
                    leaf{ind,5}, leaf{ind,6}, tol);
            end
            
            % H-Matrix distribution
            for j = 1:Nlabs
                Mh = hmxDistribution(Mh,tmp{j},ind{j},log2(Nlabs),1);
            end
            
            % Informations
            disp([' ~~> Step ',num2str(i)', ' - Elapsed time is ', ...
                num2str(hmxClock()-tps),' seconds.'])
            
            % Clean composite
            clear tmp;            
        end
    
    % Problem
    else
        error('femHmatrixParallel.m : unavailable case.');
    end
end
end


function [Mh,num] = hmxDistribution(Mh,Mp,ind,stp,num)
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
            error('femHmatrixParallel.m : unavailable case.');
        end
    end
    num = num + 1;
end
end


function tps = hmxClock()
tps = clock;
tps = tps(4)*3600 + tps(5)*60 + tps(6);
end
