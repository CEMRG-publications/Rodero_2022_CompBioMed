function scar_thinning(path,uvc,out)
clc;

%%%% Usage: 
% path: full path/name of mesh
% uvc: full path to UVC but remember there must be concordance
% out: full path/name of mesh
%

disp('Reading UVC ...')
rho=dlmread([uvc,'/COORDS_RHO.dat'],' ',0,0);
z=dlmread([uvc,'/COORDS_Z.dat'],' ',0,0);
v=dlmread([uvc,'/COORDS_V.dat'],' ',0,0);

disp('Reading points ...')
pts=dlmread([path,'.pts'],' ',1,0);

disp('Thickness calculations ...')

% Consider left ventricle
i_lv = find(v==-1);

% Consider endo and epi indices
i_epi = find(rho(i_lv)==1);
i_endo = find(rho(i_lv)==0);


h = -1*ones(size(pts,1),1);

for i=1:size(i_epi,1)
    v1 = pts(i_lv(i_epi(i)),:);
    v2 = pts(i_lv(i_endo),:);
    x = v1 - v2;
    dist = sqrt(sum(x.^2,2)); 
    h(i_lv(i_epi(i)),1) = min(dist);
 
end

disp('Write thickness matrix ...');
dlmwrite([out,'.dat'],h','delimiter',' ');

end