function Ml = hmxTimes(Ml,Mr)
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
%|    #    |   FILE       : hmxTimes.m                                    |
%|    #    |   VERSION    : 0.61                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Broadcast product with H-Matrix               |
%|  `---'  |                                                              |
%+========================================================================+

%%% H-Matrix .* H-Matrix
if isa(Ml,'hmx') && isa(Mr,'hmx')
    Ml = [];
    warning('hmxTimes.m : H-Matrix terms multiplication not yet implemented')
   

%%% Matrix .* H-Matrix -> H-Matrix
elseif isa(Mr,'hmx') 
    % H-Matrix (recursion)
    if (Mr.typ == 0)
        for i = 1:4
            if (numel(Ml)==1)
                Mr.chd{i} = hmxTimes(Ml,Mr.chd{i});
            elseif (isvector(Ml))
                Mr.chd{i} = hmxTimes(Ml(Mr.row{i}),Mr.chd{i});
            else
                error('hmxTimes.m : unavailable case')
            end
        end
        
    % Compressed leaf
    elseif (Mr.typ == 1)
        Mr.dat{1} = Ml .* Mr.dat{1};
        
    % Full leaf
    elseif (Mr.typ == 2)
        Mr.dat = Ml .* Mr.dat;

    % Unknown type
    else
        error('hmxTimes.m : unavailable case')
    end
    
    % Change H-matrix
    Ml = Mr;

    
%%% H-Matrix .* Matrix -> H-Matrix
elseif isa(Ml,'hmx') 
    % H-Matrix (recursion)
    if (Ml.typ == 0)
        for i = 1:4
            if (numel(Mr)==1)
                Ml.chd{i} = hmxTimes(Ml.chd{i},Mr);
            elseif (isvector(Mr))
                Ml.chd{i} = hmxTimes(Ml.chd{i},Mr(Ml.col{i}));            
            else
                error('hmxTimes.m : unavailable case')
            end
        end
        
    % Compressed leaf
    elseif (Ml.typ == 1)
        Ml.dat{2} = Ml.dat{2} .* Mr;
        
    % Full leaf
    elseif (Ml.typ == 2)
        Ml.dat = Ml.dat .* Mr;

    % Unknown type
    else
        error('hmxTimes.m : unavailable case')
    end
end
end
