function [A,B,flag] = ffmACA(X,Y,green,k,tol)
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
%|    #    |   FILE       : ffmACA.m                                      |
%|    #    |   VERSION    : 0.6                                           |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2019                                    |
%| ( === ) |   SYNOPSIS   : Adaptive Cross Approximation for compression  |
%|  `---'  |                of full tranfert matrix                       |
%+========================================================================+

% Partial pivoting with handle function    
row = @(i) ffmGreenKernel(X(i,:),Y,green,k).';
col = @(i) ffmGreenKernel(X,Y(i,:),green,k);

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
      
    % Compression failed
    if (n*(size(A,1)+size(B,1)) > size(A,1)*size(B,1))
        A    = [];
        B    = [];
        flag = 0;
        return
    end    
end

% B transposition lead to  A * B
B    = B.';
flag = 1;
end
