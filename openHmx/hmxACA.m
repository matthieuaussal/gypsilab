function [A,B,flag] = hmxACA(varargin)
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
%|    #    |   FILE       : hmxACA.m                                      |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Adaptative Cross Approximation, partial       |
%|  `---'  |                & total pivoting                              |
%+========================================================================+

% Total pivoting with full matrix
if (nargin == 2)
    M     = varargin{1};
    tol   = varargin{2};
    rkMax = 1e6;
    row   = @(i) full(M(i,:));
    col   = @(i) full(M(:,i));
    
% Total pivoting with full matrix and maximum rank
elseif (nargin == 3)
    M     = varargin{1};
    tol   = varargin{2};
    rkMax = varargin{3};
    row   = @(i) full(M(i,:));
    col   = @(i) full(M(:,i));
    
% Partial pivoting with handle function    
elseif (nargin == 4) 
    X     = varargin{1};
    Y     = varargin{2};
    green = varargin{3};
    tol   = varargin{4};
    rkMax = 1e6;
    row = @(i) green(X(i,:),Y).';
    col = @(i) green(X,Y(i,:));
    
% Partial pivoting with handle function and maximum rank   
elseif (nargin == 5) 
    X     = varargin{1};
    Y     = varargin{2};
    green = varargin{3};
    tol   = varargin{4};
    rkMax = varargin{5};
    row = @(i) green(X(i,:),Y).';
    col = @(i) green(X,Y(i,:));
else
    error('hmxACA : unavailable case')
end

% First row
B = row(1).';
    
% Maximum for pivoting
[~,j] = max(abs(B));
delta = B(j);
if (abs(delta) <= 1e-12)
    delta = 1e12;
end

% First column with pivot
A = col(j) ./ delta;

% Frobenius norm for the initial tensor product
An2 = (A'*A);
Bn2 = (B'*B);
Rn2 = An2 * Bn2;

% Left indices for row and columns pivot
Ir = (2:size(A,1))';
Ic = [1:j-1,j+1:size(B,1)]';

% Iterative construction
errFr = 1;
n     = 1;
while (errFr > tol)
    % Reinitialize pivot
    delta = 0;

    % Find non zeros pivot
    while abs(delta) <= 1e-12
        % No more pivots
        if isempty(Ir)
            A    = [];
            B    = [];
            flag = 0;
%             warning('hmxACA.m : no more pivots')
            return
        end
        
        % Row index
        [~,i] = max(abs(A(Ir,n)));
        
        % Compute new row
        new = row(Ir(i)).' - B*A(Ir(i),1:n).';
        
        % Column index
        [~,j] = max(abs(new(Ic)));
        
        % Pivot
        delta = new(Ic(j));
        
        % Same with minimum
        [~,iMin] = min(abs(A(Ir,n)));
        rowMin   = row(Ir(iMin)).' - B*A(Ir(iMin),1:n).';
        [~,jMin] = max(abs(rowMin(Ic)));
        deltaMin = rowMin(Ic(jMin));

        % Choose between min and max
        if abs(deltaMin) > abs(delta)
            i     = iMin;
            j     = jMin;
            new   = rowMin;
            delta = deltaMin;
        end

        % Update row indices
        Ir = [Ir(1:i-1);Ir(i+1:end)];
    end

    % Update row
    B(:,n+1) = new;
    
    % Update column
    A(:,n+1) = (col(Ic(j)) - A*B(Ic(j),1:n).') ./ delta;
    
    % Update column indices
    Ic = [Ic(1:j-1);Ic(j+1:end)];
    
    % Incrementation
    n = n + 1;

    % Recursive frobenius
    An2 = A(:,n)'*A(:,n);
    Bn2 = B(:,n)'*B(:,n);

    u  = B(:,n)'*B(:,1:n-1);
    v  = A(:,n)'*A(:,1:n-1);
    AB = v*u.';

    u  = B(:,1:n-1)'*B(:,n);
    v  = A(:,1:n-1)'*A(:,n);
    BA = u.'*v;

    Rn2 = Rn2 + AB + BA + An2*Bn2;

    % Relative Frobenius error with the residue
    errFr = sqrt(An2)*sqrt(Bn2)/sqrt(Rn2);
%     norm(A*B.' - A(:,1:n-1)*B(:,1:n-1).','fro')/norm(A*B.','fro')
%     norm(A(:,n)*B(:,n).','fro')/norm(A*B.','fro')
      
    % Compression failed
    if ( n*(size(A,1)+size(B,1)) > size(A,1)*size(B,1) ) || (n>=rkMax) 
        A    = [];
        B    = [];
        flag = 0;
%         warning('hmxACA.m : compression failed')
        return
    end    
end

% B transposition lead to  A * B
B    = B.';
flag = 1;
end
