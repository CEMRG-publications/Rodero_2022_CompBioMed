%
%HEART 4 CHAMBERS
%
%Typical call(s):
%  
%
%Variable Input:
%
%Output:
%
%Example:
%
%(c) Anton Prassl 04-2015
 
 
function mH4CsetupECGLEADS(varargin)

switch nargin
  case 0
  otherwise
    clc
    help( mfilename )
    return
end

% ---COMMENTS-------------------------------------------------------------------
% frontal view in Meshalyzer (H4C)
% z <-----0 -x
%         |
%         |
%         |
%        \/ y
%
% top view in Meshalyzer (H4C)
%         ^ x
%         |
%         |
%         |
% z <----(x) -y
%
%
% frontal view in Meshalyzer (SWG)
%         ^ z
%         |
%         |
%         |
%     -y (x) -----> x
%
% top view in Meshalyzer (SWG)
%         ^ y
%         |
%         |
%         |
%       z 0 -----> x
%
% ---DECLARATION OF CONSTANTS---------------------------------------------------
t0               = clock;
ProjPath         = sprintf('../../mesh/v1/CARP/SWGsm500um.v1.final');
carpFile         = sprintf('%s/SWGsm500um.v1', ProjPath);
ECGLeadsPath     = sprintf('%s/ECG',ProjPath);
E_Leads          = []; % 
SphereScalFactor = 1;
% ------------------------------------------------------------------------------


% IMPORT SURFACE
surf      = load_surf(carpFile);
surfNodes = unique(sort(surf.tris(:)));
XYZ       = load_pts(carpFile);


% ---PLACE EINTHOVEN LEADS------------------------------------------------------
bbox        = get_bbox(XYZ(surfNodes,:));
ctrHeart(1) = bbox.xcenter;
ctrHeart(2) = bbox.ycenter;
ctrHeart(3) = bbox.zcenter;

SphereRadius = mag(XYZ(surfNodes,[1 3]), [ctrHeart(1) ctrHeart(3)], 'sub');
SphereRadius = max(SphereRadius);

% minimum possible values to fit the Einthoven triangle
h = 3*SphereRadius; % for full coverage: 1/3h = sphere radius !
a = 2*h/sqrt(3);

if SphereScalFactor > 1
  SphereRadius = SphereRadius*SphereScalFactor;
  h            = SphereRadius*3; % for full coverage: 1/3h = sphere radius !
  a            = 2*h/sqrt(3);
end

E_Leads(end+1,:) = [ctrHeart(1)-a/2 ctrHeart(2) ctrHeart(3)+1/3*h];
E_Leads(end+1,:) = [ctrHeart(1)+a/2 ctrHeart(2) ctrHeart(3)+1/3*h];
E_Leads(end+1,:) = [ctrHeart(1)     ctrHeart(2) ctrHeart(3)-2/3*h];

mkdir(ECGLeadsPath);
save_pts(E_Leads, [ECGLeadsPath '/ecg'])
% ------------------------------------------------------------------------------


% ---PLACE WILSON LEADS---------------------------------------------------------
bbox  = get_bbox(XYZ(surfNodes,:));
w1    = ctrHeart;
w2    = ctrHeart;
w3    = ctrHeart;
w4    = ctrHeart;
w5    = ctrHeart;
w6    = ctrHeart;
w1    = w1 + [-sind( 30)*SphereRadius  -cosd( 30)*SphereRadius  0];
w2    = w2 + [ sind(  0)*SphereRadius  -cosd(  0)*SphereRadius  0];
w3    = w3 + [ sind( 30)*SphereRadius  -cosd( 30)*SphereRadius  0];
w4    = w4 + [ sind( 60)*SphereRadius  -cosd( 60)*SphereRadius  0];
w5    = w5 + [ sind( 90)*SphereRadius  -cosd( 90)*SphereRadius  0];
w6    = w6 + [ sind(120)*SphereRadius  -cosd(120)*SphereRadius  0];
% w1    = w1 + [-cosd( 60)*SphereRadius   sind( 30)*SphereRadius  0];
% w2    = w2 + [-cosd(  0)*SphereRadius  -sind(  0)*SphereRadius  0];
% w3    = w3 + [-cosd( 30)*SphereRadius  -sind( 30)*SphereRadius  0];
% w4    = w4 + [-cosd( 60)*SphereRadius  -sind( 60)*SphereRadius  0];
% w5    = w5 + [-cosd( 90)*SphereRadius  -sind( 90)*SphereRadius  0];
% w6    = w6 + [-cosd(120)*SphereRadius  -sind(120)*SphereRadius  0];
E_Leads(end+1,:) = w1;
E_Leads(end+1,:) = w2;
E_Leads(end+1,:) = w3;
E_Leads(end+1,:) = w4;
E_Leads(end+1,:) = w5;
E_Leads(end+1,:) = w6;

save_pts(E_Leads, [ECGLeadsPath '/ecg'])
% ------------------------------------------------------------------------------



% % ===TEST PLOT AREA=============================================================
% stride = 200;
% plot3(XYZ(surfNodes(1:stride:end),1),...
%       XYZ(surfNodes(1:stride:end),2),...
%       XYZ(surfNodes(1:stride:end),3), 'y')
% axis equal
% zlabel('zzz')
% ylabel('yyy')
% xlabel('xxx')
% hold on
% plot3(E_Leads(1,1), E_Leads(1,2), E_Leads(1,3), 'ro')
% line([E_Leads(1:2,1)], [E_Leads(1:2,2)], [E_Leads(1:2,3)], 'r')
% plot3(E_Leads(2,1), E_Leads(2,2), E_Leads(2,3), 'go')
% plot3(E_Leads(3,1), E_Leads(3,2), E_Leads(3,3), 'bo')
% axis tight
% % ==============================================================================




% ---EXPORT ECG LEADS AS MESHALYZER's AUXILLIARY GRID---------------------------
XYZ_leads = load_pts( [ECGLeadsPath '/ecg'] );

fid = fopen( sprintf('%s.pts_t', [ECGLeadsPath '/ecg']), 'w');
fprintf(fid, '%d\n', 3); % ? time steps

% Export nodes for the Einthoven leads as time step i
fprintf(fid, '3\n'); % #nodes for time step 1
fprintf(fid, '%.10g %.10g %.10g\n', XYZ_leads(1:3,:)'); % node

fprintf(fid, '6\n'); % #nodes for time step 2
fprintf(fid, '%.10g %.10g %.10g\n', XYZ_leads(4:9,:)'); % node

fprintf(fid, '9\n'); % #nodes for time step 3
fprintf(fid, '%.10g %.10g %.10g\n', XYZ_leads'); % node
fclose(fid);



fid = fopen( sprintf('%s.elem_t', [ECGLeadsPath '/ecg']), 'w');
fprintf(fid, '%d\n', 3); % ? time steps

% Export line elements defining Einthoven leads at time step i
fprintf(fid, '3\n');                      % #elements for time step i
fprintf(fid, 'Ln 1 0\n');                 % element
fprintf(fid, 'Ln 1 2\n');                 % element
fprintf(fid, 'Ln 0 2\n');                 % element

% Export line elements defining Einthoven leads at time step i
fprintf(fid, '5\n');                      % #elements for time step i
fprintf(fid, 'Ln 1 0\n');                 % element
fprintf(fid, 'Ln 1 2\n');                 % element
fprintf(fid, 'Ln 2 3\n');                 % element
fprintf(fid, 'Ln 3 4\n');                 % element
fprintf(fid, 'Ln 4 5\n');                 % element

% Export line elements defining Einthoven leads at time step i
fprintf(fid, '8\n');                      % #elements for time step i
fprintf(fid, 'Ln 1 0\n');                 % element
fprintf(fid, 'Ln 1 2\n');                 % element
fprintf(fid, 'Ln 2 0\n');                 % element
fprintf(fid, 'Ln 3 4\n');                 % element
fprintf(fid, 'Ln 4 5\n');                 % element
fprintf(fid, 'Ln 5 6\n');                 % element
fprintf(fid, 'Ln 6 7\n');                 % element
fprintf(fid, 'Ln 7 8\n');                 % element
fclose(fid);



data = 1;
fid  = fopen( sprintf('%s.dat_t', [ECGLeadsPath '/ecg']), 'w');
fprintf(fid, '%d\n', 3); % ? time steps

% Export line elements defining Einthoven leads at time step i
fprintf(fid, '3\n');                      % #data for time step i
fprintf(fid, '%d\n', data);               % data
fprintf(fid, '%d\n', data);               % data
fprintf(fid, '%d\n', data);               % data

% Export line elements defining Einthoven leads at time step i
fprintf(fid, '6\n');                      % #data for time step i
fprintf(fid, '%d\n', data);               % data
fprintf(fid, '%d\n', data);               % data
fprintf(fid, '%d\n', data);               % data
fprintf(fid, '%d\n', data);               % data
fprintf(fid, '%d\n', data);               % data
fprintf(fid, '%d\n', data);               % data

% Export line elements defining Einthoven leads at time step i
fprintf(fid, '9\n');                      % #data for time step i
fprintf(fid, '%d\n', data);               % data
fprintf(fid, '%d\n', data);               % data
fprintf(fid, '%d\n', data);               % data
fprintf(fid, '%d\n', data);               % data
fprintf(fid, '%d\n', data);               % data
fprintf(fid, '%d\n', data);               % data
fprintf(fid, '%d\n', data);               % data
fprintf(fid, '%d\n', data);               % data
fprintf(fid, '%d\n', data);               % data
fclose(fid);
% ------------------------------------------------------------------------------


fprintf('%c %s: %4.0f sec\n', 37, upper(mfilename), etime(clock,t0))
return