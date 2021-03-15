function mesh = Deformation_AHA(LandMarks,mesh)%load medical image ,..

addpath(genpath('./MedicalImageProcessingToolbox'))
addpath(genpath('./matlabTool'))

%We are going to assume that the 1st point is the apex, next three points are spread out of MV  the last two points
%are on the RV cusps

APEX = LandMarks.GetVertexCoordinates(1);
MV(1,:) = LandMarks.GetVertexCoordinates(2);
MV(2,:)= LandMarks.GetVertexCoordinates(3);
MV(3,:) = LandMarks.GetVertexCoordinates(4);
RV(1,:) = LandMarks.GetVertexCoordinates(5);
RV(2,:) = LandMarks.GetVertexCoordinates(6);

%in the first instance zero all points relative to apex
for i = 1:3
    MV(i,:) = MV(i,:) - APEX;
end
for i = 1:2
    RV(i,:) =RV(i,:)- APEX;
end

%Calcaulte a circle through the mitral valve points
[center,rad,v1,v2] = circlefit3d(MV(1,:),MV(2,:),MV(3,:));

%Define Rotation matrix
R = rotation2unity(center);

%Rotate points to new frame
for i = 1:2
    RV(i,:) = R*RV(i,:)';
end

theta_z = atan(-RV(2,2)/RV(2,1));
R_z = [cos(theta_z)  -sin(theta_z) 0 ; sin(theta_z) cos(theta_z) 0 ; 0 0 1];

for i = 1:2
    RV(i,:) = R_z*RV(i,:)';
end
R = R_z*R;

for i = 1:3
    MV(i,:) = R*MV(i,:)';
end

[center,~,~] = circlefit3d(MV(1,:),MV(2,:),MV(3,:));

%Rotate mesh into new coordiantes
for i = 1:length(mesh.points)
    mesh.points(i,:) = (R*(mesh.points(i,:) - APEX)')';
end

%Creat label indexes
% oLab = [ 2 3 4 5 6 1;
%          8 9 10 11 12 7;
%          14 15 16 13 0 0;
%          17 0 0 0 0 0];
oLab = [ 5 4 3 2 1 6;
    11 10 9 8 7 12;
    15 14 13 16 0 0;
    17 0 0 0 0 0];

apical_segs_rotate = pi/4;

if (center(3) < 0)
    mesh.points(:,3) = -mesh.points(:,3);
    oLab = [ 2 3 4 5 6 1;
             8 9 10 11 12 7;
             14 15 16 13 0 0;
             17 0 0 0 0 0];
%     apical_segs_rotate = -apical_segs_rotate;
end

%Create AHA segments from new coordinates
%Find z range
Zmax = max(mesh.points(:,3));
Zmin = 0;
RangeZ = (Zmax - Zmin)*0.8;

TOP= RangeZ+Zmin ;
BASE = RangeZ*2/3+Zmin;
MID =  RangeZ*1/3+Zmin;

%Angle RV cusp 2
RVangle1 = abs(atan2(RV(1,2),RV(1,1)));
RVangle2 = abs(atan2(RV(2,2),RV(2,1)));

%assuming RVangle1< RV angle 2
sepA = abs(RVangle2 - RVangle1)/2;
freeA = (2*pi- abs(RVangle2 - RVangle1))/4;

%point angles
angle = (atan2(mesh.points(:,2),mesh.points(:,1)))-RVangle1;
angle = (angle).*(angle>0) + (2*pi+angle).*(angle<0);

%Base points
Bindex = (mesh.points(:,3)>=BASE).* (mesh.points(:,3)<TOP);
Mindex = (mesh.points(:,3)>=MID).* (mesh.points(:,3)<BASE);
Aindex = (mesh.points(:,3)>=Zmin).* (mesh.points(:,3)<MID);
AAindex = (mesh.points(:,3)<Zmin);

LABEL = zeros(1,mesh.npoints);

%create scaler 
SCAL = 1;

%loop over three layers
for k = 1:4
    if(k==1)
        INDEX = Bindex;
        nSec = 6;
        WID = [freeA freeA sepA sepA freeA freeA];
        if (center(3) < 0)
             WID = [sepA sepA freeA freeA freeA freeA];
        end
        Csec = 0;
    elseif (k==2)
        INDEX = Mindex;
        nSec = 6;
        WID = [freeA freeA sepA sepA freeA freeA];
        if (center(3) < 0)
             WID = [sepA sepA freeA freeA freeA freeA];
        end
        Csec = 0;
    elseif(k==3)
        INDEX = Aindex;
        nSec = 4;
        WID = [pi/2 pi/2 pi/2 pi/2 ];
        Csec = sepA - apical_segs_rotate;
    elseif(k==4)
        INDEX = AAindex;
        nSec = 1;
        WID = [2*pi];
        Csec = 0;
    end
    
    for i = 1:nSec
        Upper = Csec +WID(i);
        Lower = Csec;
        
        LABEL(logical((angle<Upper).*(angle>=Lower).*INDEX)) = oLab(k,i)*SCAL;
        if Lower <0
            LABEL(logical((angle<2*pi).*(angle>=(2*pi + Lower)).*INDEX)) = oLab(k,i)*SCAL;
        end
        if Upper > 2*pi
            LABEL(logical((angle<Upper -2*pi).*(angle>=0).*INDEX)) = oLab(k,i)*SCAL;
        end
        Csec = Csec + WID(i);
    end
end

%assign labels to mesh
mesh.attributes(1).attribute_array = LABEL;
mesh.attributes(1).nelements = length(LABEL);
mesh.addVertexAttribute(LABEL,'AHA_Label');



