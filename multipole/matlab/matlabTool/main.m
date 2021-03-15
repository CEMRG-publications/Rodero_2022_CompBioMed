close all
clear all
addpath(genpath('/home/or15/Install/MedicalImageProcessingToolbox'))

LM_FILE = '../LMM.vtk';

lm = read_vtkMesh(LM_FILE);

%load and label reference mesh
Rmesh =  read_vtkMesh(sprintf('../../proc_CT/segmentation_ascii-%i.vtk',0));
Dmesh = Deformation_AHA(lm,Rmesh);
Dmesh.calcualteArea;
for k = 1:16
    Rarea(k) = sum(Dmesh.Area(Dmesh.attributes(2).attribute_array==k));
end
    

%Show mesh with AHA segmentseg_ascii
viewMesh(Dmesh,'labelColor',1)

%We want to load in all 10 meshes and then calcaulte the area
for i = 1:10;
    MeshStack(i) = read_vtkMesh(sprintf('../../proc_CT/segmentation_ascii-%i.vtk',i-1));
    MeshStack(i).calcualteArea;
    SQUEEZE(i,:) = (MeshStack(i).Area -  Dmesh.Area)./Dmesh.Area;
    for k = 1:16
        Rsqueez(i,k) = sum(MeshStack(i).Area(Dmesh.attributes(2).attribute_array==k).*SQUEEZE(i,Dmesh.attributes(2).attribute_array==k))/Rarea(k);
    end
   
end

%Calcualte area wieghted change in each region

figure
hold on
plot(0:10,[Rsqueez(:,1:6) ; Rsqueez(1,1:6)],'--')
plot(0:10,[Rsqueez(:,7:12) ; Rsqueez(1,7:12)],'.-')
plot(0:10,[Rsqueez(:,13:16) ; Rsqueez(1,13:16)],'-')



legend('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','Location','NorthWest')