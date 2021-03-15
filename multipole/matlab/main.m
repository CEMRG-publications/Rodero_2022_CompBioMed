clear
addpath(genpath('./MedicalImageProcessingToolbox'))
addpath(genpath('./matlabTool'))

folder = '/media/crg17/Seagate Backup Plus Drive/CT_cases/h_case20/meshing/1000um/BiV';

LM_FILE = fullfile(folder,'LM.vtk');

lm = read_vtkMesh(LM_FILE);

%load and label reference mesh
Rmesh =  read_vtkMesh(fullfile(folder,'BiV.LV.vtk'));
Dmesh = Deformation_AHA_new(lm,Rmesh);
% write_vtkMesh(fullfile(folder,'AHA.vtk'), Dmesh);

% Write labels separately in a .dat file to map it on the elements
FileOutID = fullfile(folder,'AHA.dat');
FileOut = fopen(FileOutID,'wt');
for i = 1:length(Dmesh.attributes(1).attribute_array)
    fprintf(FileOut,'%d\n',Dmesh.attributes(1).attribute_array(i));
end
fclose(FileOut);
% 
% load handel.mat;
% sound(y(1:2e4), Fs);

% command=['/home/common/bin/GlVTKConvert -m ', fullfile(folder,'BiV.LV'),...
%                      ' -n ', fullfile(folder,'AHA.dat'),...
%                      ' -o ', fullfile(folder,'AHA_02')];
% system(command)

% generate default files to remap the AHA segments onto the four-chamber
% geometry
% nPts = 462282;
% FileOutID = fullfile(folder,'AHA_fourChamber.dat');
% FileOut = fopen(FileOutID,'wt');
% for i = 1:nPts
%     fprintf(FileOut,'0\n');
% end
% fclose(FileOut);

