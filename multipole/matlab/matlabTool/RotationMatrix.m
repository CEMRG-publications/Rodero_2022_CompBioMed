function R = RotationMatrix(u,theta)

R(1,1) = u(1)^2 +cos(theta) * (1 - u(1)^2);
R(1,2) = (1 -cos(theta)) * u(1) * u(2) - u(3) *sin(theta);
R(1,3) = (1 -cos(theta)) * u(1) * u(3) + u(2) *sin(theta);

R(2,1) = (1 -cos(theta)) * u(1) * u(2) + u(3) *sin(theta);
R(2,2) = u(2)^2 +cos(theta) * (1 - u(2)^2);
R(2,3) = ( 1 - cos(theta)) * u(2) * u(3) - u(1) *sin(theta);

R(3,1) = ( 1 - cos(theta)) * u(1) * u(3) - u(2) *sin(theta);
R(3,2) = ( 1 - cos(theta)) * u(2) * u(3) + u(1) *sin(theta);
R(3,3) = u(3)^2 +cos(theta) * (1 - u(3)^2);

end