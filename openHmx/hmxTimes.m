function Mh = hmxTimes(Ml,Mr)
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
%|    #    |   FILE       : hmxTimes.m                                    |
%|    #    |   VERSION    : 0.32                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 25.12.2017                                    |
%| ( === ) |   SYNOPSIS   : Scalar product with H-Matrix                  |
%|  `---'  |                                                              |
%+========================================================================+

%%% H-Matrix .* H-Matrix
if isa(Ml,'hmx') && isa(Mr,'hmx')
    Mh = [];
    warning('hmxTimes.m : H-Matrix terms multiplication not yet implemented')
   
    
%%%  H-Matrix .* scal   ||   scal .* H-Matrix -> H-Matrix
else
    % Input analysis
    if isa(Ml,'hmx')
        Mh = Ml;
        x  = Mr;
    elseif isa(Mr,'hmx')
        Mh = Mr;
        x  = Ml;
    else
        error('hmxTimes.m : unavailable case')
    end
    
    % H-Matrix (recursion)
    if (Mh.typ == 0)
        for i = 1:4
            Mh.chd{i} = hmxTimes(x,Mh.chd{i});
        end
        
    % Compressed leaf
    elseif (Mh.typ == 1)
        Mh.dat{1} = x .* Mh.dat{1};
        
    % Full leaf
    elseif (Mh.typ == 2)
        Mh.dat = x .* Mh.dat;

    % Unknown type
    else
        error('hmxTimes.m : unavailable case')
    end
end
end
