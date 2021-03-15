function mesh = Deformation_AHA_new(LandMarks,mesh)%load medical image ,..

addpath(genpath('./MedicalImageProcessingToolbox'))
addpath(genpath('./matlabTool'))

warning('The first RV point must be anterior');

%We are going to assume that the 1st point is the apex, next three points are spread out of MV  the last two points
%are on the RV cusps

APEX = LandMarks.GetVertexCoordinates(1);
MV(1,:) = LandMarks.GetVertexCoordinates(2);
MV(2,:)= LandMarks.GetVertexCoordinates(3);
MV(3,:) = LandMarks.GetVertexCoordinates(4);
RV(1,:) = LandMarks.GetVertexCoordinates(5);
RV(2,:) = LandMarks.GetVertexCoordinates(6);

% Define plane 
n = cross((MV(2,:)-MV(1,:))./norm(MV(2,:)-MV(1,:)),(MV(3,:)-MV(1,:))./norm(MV(3,:)-MV(1,:)));
direction = dot(APEX-MV(1,:),n);

if (direction >=0)
    n = -n./norm(n);
else
    n = n./norm(n);
end

% Define points along the height
c = mean(MV,1);
long_axis = norm(c-APEX);
n_APEX = (c-APEX)./long_axis;

c_midbase = c - n_APEX*long_axis/3;
c_mid = c - n_APEX*long_axis*2/3;

LABEL = zeros(size(mesh.points(:,1)));
for i = 1:length(mesh.points)
    if(dot(mesh.points(i,:)-c,n)>0)
        LABEL(i) = 0;
    elseif((dot(mesh.points(i,:)-c_midbase,n)>0) && (dot(mesh.points(i,:)-c,n)<=0))
        LABEL(i) = 1;
    elseif((dot(mesh.points(i,:)-c_mid,n)>0) && (dot(mesh.points(i,:)-c_midbase,n)<=0))  
        LABEL(i) = 7;
    elseif((dot(mesh.points(i,:)-APEX,n)>0) && (dot(mesh.points(i,:)-c_mid,n)<=0))    
        LABEL(i) = 13;
    elseif(dot(mesh.points(i,:)-APEX,n)<=0)
        LABEL(i) = 17;
    end
end

% COG_apex = c - n_APEX*long_axis*5/6;
COG_apex = c - n_APEX*long_axis*2/3;

%% planes for rv junctions
rvj1_base = RV(1,:) - dot(RV(1,:) - c,n)*n;
rvj2_base = RV(2,:) - dot(RV(2,:) - c,n)*n;
midSept_base = 0.5*(rvj1_base+rvj2_base);

rvj1_mid = rvj1_base - norm(c - c_midbase)*n_APEX;
rvj2_mid = rvj2_base - norm(c - c_midbase)*n_APEX;

midSept_mid = midSept_base - norm(c - c_midbase)*n_APEX;

n_rvj1 = cross((c-c_midbase)./norm(c-c_midbase),(c_midbase-rvj1_mid)./norm(c_midbase-rvj1_mid)); 
n_rvj2 = cross((c-c_midbase)./norm(c-c_midbase),(c_midbase-rvj2_mid)./norm(c_midbase-rvj2_mid)); 
n_midSept = cross((c-c_midbase)./norm(c-c_midbase),(c_midbase-midSept_mid)./norm(c_midbase-midSept_mid)); 

%% planes for lv 
rvj1_angle = acos(dot((rvj1_base-c)./norm(rvj1_base-c),(rvj2_base-c)./norm(rvj2_base-c)));
lv_angle = (2*pi - rvj1_angle)/4;

%rotation axis
long_axis_v = (c-APEX)./norm(c-APEX);
R = RotationMatrix(long_axis_v,-lv_angle);

p45_midbase_v = R*normr(rvj2_mid-c_midbase)';
p45_midbase = c_midbase + p45_midbase_v' * norm(rvj2_mid-c_midbase);
n_45 = cross((c-c_midbase)./norm(c-c_midbase),(c_midbase-p45_midbase)./norm(c_midbase-p45_midbase)); 
Deformation_AHA.m
p56_midbase_v = R*normr(p45_midbase-c_midbase)';
p56_midbase = c_midbase + p56_midbase_v' * norm(p45_midbase-c_midbase);
n_56 = cross((c-c_midbase)./norm(c-c_midbase),(c_midbase-p56_midbase)./norm(c_midbase-p56_midbase)); 

p61_midbase_v = R*normr(p56_midbase-c)';
p61_midbase = c_midbase + p61_midbase_v' * norm(p56_midbase-c_midbase);
n_61 = cross((c-c_midbase)./norm(c-c_midbase),(c_midbase-p61_midbase)./norm(c_midbase-p61_midbase));  

R_adjust = RotationMatrix(long_axis_v,pi/4);
midSept_apex = midSept_mid - norm(c_midbase - COG_apex)*n_APEX;
p1314_v = R_adjust*normr(midSept_apex-COG_apex)';
p1314_apex = COG_apex + p1314_v' * norm(midSept_apex-COG_apex);
n_1314 = cross((c_midbase-COG_apex)./norm(c_midbase-COG_apex),(COG_apex-p1314_apex)./norm(COG_apex-p1314_apex));

R_apex = RotationMatrix(long_axis_v,-pi/2);

p1415_v = R_apex*normr(p1314_apex-COG_apex)';
p1415_apex = COG_apex + p1415_v' * norm(p1314_apex-COG_apex);
n_1415 = cross((c_midbase-COG_apex)./norm(c_midbase-COG_apex),(COG_apex-p1415_apex)./norm(COG_apex-p1415_apex)); 

p1516_v = R_apex*normr(p1415_apex-COG_apex)';
p1516_apex = COG_apex + p1516_v' * norm(p1415_apex-COG_apex);
n_1516 = cross((c_midbase-COG_apex)./norm(c_midbase-COG_apex),(COG_apex-p1516_apex)./norm(COG_apex-p1516_apex)); 

p1713_v = R_apex*normr(p1516_apex-COG_apex)';
p1713_apex = COG_apex + p1713_v' * norm(p1516_apex-COG_apex);
n_1713 = cross((c_midbase-COG_apex)./norm(c_midbase-COG_apex),(COG_apex-p1713_apex)./norm(COG_apex-p1713_apex)); 

for i = 1:length(mesh.points)
    if((LABEL(i)>=1) && (LABEL(i)<=7))
        if((dot(mesh.points(i,:)-rvj1_mid,n_rvj1)>0) && (dot(mesh.points(i,:)-midSept_mid,n_midSept)<=0))
            LABEL(i) = LABEL(i) + 1;
        elseif((dot(mesh.points(i,:)-midSept_mid,n_midSept)>0) && (dot(mesh.points(i,:)-rvj2_mid,n_rvj2)<=0))
            LABEL(i) = LABEL(i) + 2;
        elseif((dot(mesh.points(i,:)-rvj2_mid,n_rvj2)>0) && (dot(mesh.points(i,:)-p45_midbase,n_45)<=0))
            LABEL(i) = LABEL(i) + 3;
        elseif((dot(mesh.points(i,:)-p45_midbase,n_45)>0) && (dot(mesh.points(i,:)-p56_midbase,n_56)<=0))
            LABEL(i) = LABEL(i) + 4;
        elseif((dot(mesh.points(i,:)-p56_midbase,n_56)>0) && (dot(mesh.points(i,:)-p61_midbase,n_61)<=0))
            LABEL(i) = LABEL(i) + 5;
        end
    elseif((LABEL(i)>7) && (LABEL(i)<17))
        if((dot(mesh.points(i,:)-p1314_apex,n_1314)>0) && (dot(mesh.points(i,:)-p1415_apex,n_1415)<=0))
            LABEL(i) = LABEL(i) + 1;
        elseif((dot(mesh.points(i,:)-p1415_apex,n_1415)>0) && (dot(mesh.points(i,:)-p1516_apex,n_1516)<=0))
            LABEL(i) = LABEL(i) + 2;
        elseif((dot(mesh.points(i,:)-p1516_apex,n_1516)>0) && (dot(mesh.points(i,:)-p1713_apex,n_1713)<=0))
            LABEL(i) = LABEL(i) + 3;
        end
    end
end

% assign labels to mesh
mesh.attributes(1).attribute_array = LABEL;
mesh.attributes(1).nelements = length(LABEL);
mesh.addVertexAttribute(LABEL,'AHA_Label');

