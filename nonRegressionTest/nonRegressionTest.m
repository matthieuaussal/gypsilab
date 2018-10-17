%+========================================================================+
%|                                                                        |
%|            This script uses the GYPSILAB toolbox for Matlab            |
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
%|    #    |   FILE       : nonRegressionTest.m                           |
%|    #    |   VERSION    : 0.41                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 01.04.2018                                    |
%| ( === ) |   SYNOPSIS   : run all non regression test                   |
%|  `---'  |                                                              |
%+========================================================================+

clear all
close all
clc

% Mesh management
cd('meshManagement/')
run('nrtMshBitree.m')
run('nrtMshClean.m');
run('nrtMshCube.m');
run('nrtMshOctree.m')
run('nrtMshRead.m');
run('nrtMshRefine.m');
run('nrtMshSegment.m');
run('nrtMshSquare.m');

% Domain quadrature
cd('../domainQuadrature/')
run('nrtDom1D.m');
run('nrtDomHalfSquare.m');
run('nrtDomIntSing.m');
run('nrtDomND.m');
run('nrtDomTetra.m');
run('nrtDomTrace.m');
run('nrtDomTriangle.m');

% Finite element
cd('../finiteElement/')
run('nrtFemAllWithShuffle.m');
run('nrtFemContinuity1.m');
run('nrtFemContinuity2.m');
run('nrtFemConvergence.m');
run('nrtFemDirichletCube.m');
run('nrtFemDirichletDisk.m');
run('nrtFemDirichletSquare.m');
run('nrtFemDirichletString.m');
run('nrtFemElasticite2D')
run('nrtFemJunction.m');
run('nrtFemLaplace.m');
run('nrtFemOperators.m');
run('nrtFemRwgNed.m');
run('nrtFemWave1D.m');
run('nrtFemWave2D.m');

% Hierarchical matrices
cd('../hierarchicalMatrix')
run('nrtHmxAlgebra.m');
run('nrtHmxBuilder.m');
run('nrtHmxBuilderFem.m');
run('nrtHmxCompareLU.m');
run('nrtHmxCompressorPartial.m');
run('nrtHmxCompressorsBox.m');
run('nrtHmxCompressorTotal.m');

% Scattering
cd('../scattering')
run('nrtHelmholtzCalderon.m')
run('nrtHmxHelmholtzS.m');
run('nrtHmxHelmholtzD.m');
run('nrtHmxHelmholtzDt.m');
run('nrtHmxHelmholtzH.m');
run('nrtHmxHelmholtzSDrad.m');
run('nrtHmxHelmholtzBWdir.m');
run('nrtHmxHelmholtzBWneu.m');
run('nrtHmxMaxwellT.m');
run('nrtHmxMaxwellNxK.m');
run('nrtHmxMaxwellCFIE.m');

% Fem-Bem dielectrique
cd('../femBemDielectrique')
run('nrtHmxFemBemEFIE.m');
run('nrtHmxFemBemEFIEhalf.m');
run('nrtHmxFemBemCFIE.m');
run('nrtHmxFemBemCFIEhalf.m');

% Inverse problem
cd('../inverseProblem')
run('nrtIpbHelmholtz.m');
run('nrtIpbHelmholtz0.m');

% Ray-tracing
cd('../rayTracing');
run('nrtRayCube.m');
run('nrtRayFabryPerot.m');
run('nrtRayLabyrinthe.m');
run('nrtRaySource.m');
run('nrtRaySphere.m');
run('nrtRayTheatre');

% End
disp('~~> Non Regresion Test are done. Michto gypsilab !');


