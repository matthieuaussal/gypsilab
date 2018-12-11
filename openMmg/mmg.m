function mesh = mmg(mesh,hmin,hmax)
% Copyright (c) 20015-2018, Matthieu Aussal, Ecole Polytechnique       
% LGPL Lesser General Public License v3.0. 
% Remeshing using level-set from Mmg tools : https://www.mmgtools.org           

% Current and mmg directory
here  = pwd;
there = which('mmg.m');
there = there(1:end-6);

% Move to mmg directory
cd(there)

% Security
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

% Convert input to command
file = '-in original.msh -out refined.msh ';  % mesh temporary files 
hmin = ['-hmin  ',num2str(hmin),' '];   % edge min
hmax = ['-hmax  ',num2str(hmax),' '];   % edge max
info = ['-v ',num2str(0),' ']; % 0 (no infos) to 6 (everything), 1 default 

% Execute mmg binaries
command = [os,bin,file,hmin,hmax,info]
system(command);  

% Read refined mesh
mesh = msh('refined.msh');

% Clean meshes
delete('original.msh')
delete('refined.msh')

% Back to current directory
cd(here)
end
