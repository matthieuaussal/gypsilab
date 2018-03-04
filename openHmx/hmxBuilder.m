function Mh = hmxBuilder(X,Y,green,tol)
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
%|    #    |   FILE       : hmxBuilder.m                                  |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Particles builder with low-rank approximation |
%|  `---'  |                for full, sparse and handle function          |
%+========================================================================+

% Initialisation
Mh = hmx(X,Y,tol);

% Compression for far distances
if ~issparse(green) && hmxFar(Mh)
    % ACA with partial pivoting for handle function
    if isa(green,'function_handle')
        [A,B,flag] = hmxACA(X,Y,green,tol);
        
    % ACA with total pivoting for full matrix
    else
        [A,B,flag] = hmxACA(green,tol);
    end
    
% Compression for empty box
elseif issparse(green) && (nnz(green) == 0)
    A    = zeros(size(green,1),0);
    B    = zeros(0,size(green,2));
    flag = 1;
    
% No compression    
else
    flag = 0;
end


%%% Compression
if flag
    Mh.dat = {A,B};
    Mh.typ = 1;
    
    
%%%% Full or sparse for smallest box (stopping criterion)
elseif sum(Mh.dim < 100)
    % Handle function
    if isa(green,'function_handle')
        [I,J]  = ndgrid(1:size(X,1),1:size(Y,1));
        Mh.dat = green(X(I,:),Y(J,:));
        Mh.dat = reshape(Mh.dat,size(X,1),size(Y,1));
    
    % Full or sparse matrix
    else
        Mh.dat = green;
    end
    
    % Type
    Mh.typ = 2;
    
    
%%% H-Matrix (recursion)
else
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
        % Coordinates
        Xi = X(Mh.row{i},:);
        Yi = Y(Mh.col{i},:);
        
        % Partial pivoting
        if isa(green,'function_handle')
            Mh.chd{i} = hmxBuilder(Xi,Yi,green,tol);
            
        % Total pivoting    
        elseif isnumeric(green)
            Mi = green(Mh.row{i},Mh.col{i});
            Mh.chd{i} = hmxBuilder(Xi,Yi,Mi,tol);
            
        else
            error('hmxBuilder.m : unavailable case')
        end
    end
    
    % Type
    Mh.typ = 0; 
    
    % Fusion
    Mh = hmxFusion(Mh);    
end
end
