function mesh = mmg(varargin)
% Copyright (c) 20015-2018, Matthieu Aussal, Ecole Polytechnique       
% LGPL Lesser General Public License v3.0. 
% Remeshing using level-set from Mmg tools : https://www.mmgtools.org           

% Current and mmg directory
here  = pwd;
there = which('mmg.m');
there = there(1:end-6);

% Move to mmg directory
cd(there)

% Check mesh (3D nodes, triangle or tetra) 
mesh = varargin{1};
if (size(mesh.vtx,2) ~= 3) || (size(mesh.elt,2) < 3) || (size(mesh.elt,2) > 4)
    error('mmg.m : unavailable case')
end

% Write original mesh in .msh format
mshWriteMsh('original.msh',mesh);

% Operating system
if ismac
    os = './mac_';
elseif ispc
    os = './win_';
elseif isunix
    os = './uni_';
else
    disp('mmg.m : unavailable case')
end

% Element form
if (size(mesh.elt,2) == 3)
    bin = 'mmgs_O3';
elseif (size(mesh.elt,2) == 4)
    bin = 'mmg3d_O3';    
end

% Add extension for windows only
if ispc
    bin = [bin,'.exe '];
else
    bin = [bin,' '];
end

% Define temporary mesh files
file = '-in original.msh -out refined.msh ';  % mesh temporary files 

% Convert input to command
if (nargin==1)
    stp = mesh.stp;
    opt = ['-hsiz ',num2str(stp(3)),' '];  
elseif (nargin==2)
    opt = ['-hmin  ',num2str(varargin{2}(1)),' -hmax  ',num2str(varargin{2}(2)),' '];   
else
    error('mmg.m : unavailable case');
end
    
% Define verbose level : 0 (no infos), 1 default to 6 (everything)
info = ['-v ',num2str(0),' ']; 

% Execute mmg binaries
command = [os,bin,file,opt,info];
system(command);  

% Read refined mesh
mesh = msh('refined.msh');

% Clean meshes
delete('original.msh')
delete('refined.msh')

% Back to current directory
cd(here)
end
