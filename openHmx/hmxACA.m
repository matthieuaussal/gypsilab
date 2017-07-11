function [A,B,flag] = hmxACA(X,Y,green,tol)
%+========================================================================+
%|                                                                        |
%|               OPENHMX, H-MATRIX COMPRESSION AND ALGEBRA                |
%|              openHmx is part of GYPSYLAB toolbox - v0.20               |
%|                                                                        |
%| Copyright (c) 20015-2017, Ecole polytechnique, all rights reserved.    |
%| Licence Creative Commons BY-NC-SA 4.0, Attribution, NonCommercial and  |
%| ShareAlike (see http://creativecommons.org/licenses/by-nc-sa/4.0/).    |
%| This software is the property from Centre de Mathematiques Appliquees  |
%| de l'Ecole polytechnique, route de Saclay, 91128 Palaiseau, France.    |
%|                                                            _   _   _   |
%| Please acknowledge the GYPSILAB toolbox in programs       | | | | | |  |
%| or publications in which you use the code. For openHmx,    \ \| |/ /   |
%| we suggest as reference :                                   \ | | /    |
%| [1] : www.cmap.polytechnique.fr/~aussal/gypsilab             \   /     |
%| [2] : 13th International Conference on Mathematical           | |      |
%| and Numerical Aspects of Wave Propagation, University of      | |      |
%| Minnesota, may 2017. "OpenHmX, an open-source H-Matrix        | |      |
%| toolbox in Matlab".                                           | |      |
%|_______________________________________________________________|_|______|
%| Author(s)  : Matthieu Aussal - CMAP, Ecole polytechnique               |
%| Creation   : 14.03.17                                                  |
%| Last modif : 21.06.17                                                  |
%| Synopsis   : Adaptative Cross Approximation, partial & total pivoting  |
%+========================================================================+

% Partial or total pivoting
if isnumeric(green)
    mtx = 1;
else 
    mtx = 0;
end

% Dimensions
Nx = size(X,1);
Ny = size(Y,1);

% First row
if mtx
    B = green(1,:).';
else
    B = green(X(1,:),Y);
end

% Maximum for pivoting
[delta,j] = max(B);

% First column
if abs(delta) <= 10*eps(B(1))
    A = zeros(Nx,1,class(B));
elseif mtx
    A = green(:,1) ./ delta;
else
    A = green(X,Y(j,:)) ./ delta;
end

% Sorted interactions
X0     = mean(X,1);
Y0     = mean(Y,1);
[~,ix] = sort( (X-ones(Nx,1)*X0) * (Y0-X0)' );
[~,iy] = sort( (Y-ones(Ny,1)*Y0) * (X0-Y0)' );

% Length compatibility
if Nx > Ny
    iy = iy(mod(0:Nx-1,Ny)'+1);
else
    ix = ix(mod(0:Ny-1,Nx)'+1);
end

% Add random interactions
ix  = [ix ; ceil(Nx*rand(max(Nx,Ny),1))];
iy  = [iy ; ceil(Ny*rand(max(Nx,Ny),1))];

% Reference
if mtx
    ref = green(sub2ind(size(green),ix,iy));
else
    ref = green(X(ix,:),Y(iy,:));
end
nrf = norm(ref,'inf');

% Numerical solution
sol = A(ix) .* B(iy);

% Construction
n = 1;
while (norm(ref-sol,'inf')/nrf > tol)
    % Initialisation
    ind   = (1:Nx)';
    i     = 0;
    delta = 0;
    
    % Find non zeros pivot
    while abs(delta) <= 10*eps(B(1))
        % Extract void pivot
        if i
            ind = ind([1:i-1,i+1:end]);
            if isempty(ind)
                A = [];
                B = [];
                flag = 0;
                return
            end
        end
        
        % Row index
        [~,i] = max(A(ind,n));
        
        % Compute new row
        if mtx
            row = green(ind(i),:).' - B*A(ind(i),1:n).';
        else
            row = green(X(ind(i),:),Y) - B*A(ind(i),1:n).';
        end
        
        % Column index
        [delta,j] = max(row);
    end
    
    % Update row
    B(:,n+1) = row;
    
    % Compute new column
    if mtx
        A(:,n+1) = (green(:,j) - A*B(j,1:n).') ./ delta;
    else
        A(:,n+1) = (green(X,Y(j,:)) - A*B(j,1:n).') ./ delta;
    end
    
    % Incrementation
    n = n + 1;
    
    % Accuracy of the compression
    sol = sol + A(ix,n) .* B(iy,n);
    
    % Compression failed
    if (numel(A)+numel(B) > Nx*Ny)
        A = [];
        B = [];
        flag = 0;
        return
    end
end

% B transposition lead to  A * B
B    = B.';
flag = 1;
end
