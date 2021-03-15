function R = rotation2unity(v)

theta_x = atan(v(2)/v(3));
R_x = [1 0 0 ; 0 cos(theta_x) -sin(theta_x); 0 sin(theta_x) cos(theta_x)];
v_x = R_x*v';
theta_y = atan(-v_x(1)/v_x(3));
R_y = [cos(theta_y) 0 sin(theta_y); 0 1 0; -sin(theta_y) 0 cos(theta_y)];
R = R_y*R_x;