function mesh = Deformation_AHA(LandMarks,mesh)%load medical image ,..

addpath(genpath('/home/or15/Install/MedicalImageProcessingToolbox'))

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

%Calcaulte a circle through teh mitral valve points
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

Rcenter = R*center';
for i = 1:3
    MV(i,:) = R*MV(i,:)';
end



hold on;
for i =1:3
    plot3(MV(i,1),MV(i,2),MV(i,3),'bo');
end

for i =1:2
    plot3(RV(i,1),RV(i,2),RV(i,3),'ro');
end

plot3(0,0,0,'go');

plot3([0,Rcenter(1)],[0,Rcenter(2)],[0,Rcenter(3)]);

plot3([0,RV(1,1)],[0,RV(1,2)],[RV(1,3),RV(1,3)]);

%Rotate mesh into new coordiantes
for i = 1:length(mesh.points)
    mesh.points(i,:) = (R*(mesh.points(i,:) - APEX)')';
end

%update centres
mesh.calcualteArea;

 viewMesh(mesh)
%Create AHA segments from new coordinates
%Find z range
Zmax = max(mesh.points(:,3));
Zmin = min(mesh.points(:,3));
RangeZ = (Zmax - Zmin)*0.8;

TOP= RangeZ+Zmin ;
BASE = RangeZ*2/3+Zmin;
MID =  RangeZ*1/3+Zmin;

%Angle RV cusp 2
RVangle1 = atan2(RV(1,2),RV(1,1));
RVangle2 = atan2(RV(2,2),RV(2,1));

%assuming RVangle1< RV angle 2
sepA = (RVangle2 - RVangle1)/2;
freeA = (2*pi- (RVangle2 - RVangle1))/4;

%point angles
angle = (atan2(mesh.points(:,2),mesh.points(:,1)))-RVangle1;
angle = (angle).*(angle>0) + (2*pi+angle).*(angle<0);

%Centre angles
Cangle = (atan2(mesh.centres(:,2),mesh.centres(:,1)))-RVangle1;
Cangle = (Cangle).*(Cangle>0) + (2*pi+Cangle).*(Cangle<0);

%Base points
Bindex = (mesh.points(:,3)>=BASE).* (mesh.points(:,3)<TOP);
Mindex = ((mesh.points(:,3)>=MID).* (mesh.points(:,3)<BASE) );
Aindex = (mesh.points(:,3)<MID);

cBindex = (mesh.centres(:,3)>=BASE).* (mesh.centres(:,3)<TOP);
cMindex = ((mesh.centres(:,3)>=MID).* (mesh.centres(:,3)<BASE) );
cAindex = (mesh.centres(:,3)<MID);

LABEL = zeros(1,mesh.npoints);
cLABEL = zeros(1,mesh.ntriangles);

%Creat label indexes
oLab = [ 2 3 4 5 6 1;
    8 9 10 11 12 7
    14 15 16 13 0 0];

%create scaler 
SCAL = length(angle)/16;

%loop over three layers
for k = 1:3
    if(k==1)
        INDEX = Bindex;
         cINDEX = cBindex;       
        nSec = 6;
        WID = [sepA sepA freeA freeA freeA freeA];
        Csec = 0;
    elseif (k==2)
        INDEX = Mindex;
         cINDEX = cMindex; 
        nSec = 6;
        WID = [sepA sepA freeA freeA freeA freeA];
        Csec = 0;
    elseif(k==3)
        INDEX = Aindex;
        cINDEX = cAindex; 
        nSec = 4;
        WID = [pi/2 pi/2 pi/2 pi/2 ];
        Csec = sepA - pi/4;
    end
    
    for i = 1:nSec
        Upper = Csec +WID(i);
        Lower = Csec;
        
        LABEL(logical((angle<Upper).*(angle>=Lower).*INDEX)) = oLab(k,i)*SCAL;
        cLABEL(logical((Cangle<Upper).*(Cangle>=Lower).*cINDEX)) = oLab(k,i);       
        if Lower <0
            LABEL(logical((angle<2*pi).*(angle>=(2*pi + Lower)).*INDEX)) = oLab(k,i)*SCAL;
            cLABEL(logical((Cangle<2*pi).*(Cangle>=(2*pi + Lower)).*cINDEX)) = oLab(k,i);         
        end
        if Upper > 2*pi
            LABEL(logical((angle<Upper -2*pi).*(angle>=0).*INDEX)) = oLab(k,i)*SCAL;
            cLABEL(logical((Cangle<Upper -2*pi).*(Cangle>=0).*cINDEX)) = oLab(k,i);           
        end
        Csec = Csec + WID(i);
    end
end

%assign labels to mesh
mesh.attributes(1).attribute_array = LABEL;
mesh.attributes(2).attribute_array = cLABEL;
mesh.addVertexAttribute(LABEL,'AHA_Label');



