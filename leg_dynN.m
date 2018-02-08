function [xdot, xdot_x, xdot_u] = leg_dynN(x, u)

% x is state vector of dim 6:
% theta 1
% theta 2
% theta 3
% dtheta 1
% dtheta 2
% dtheta 3

% u is control vector of dim 3:
% torque 1 (along y-axis)
% torque 2 (along x-axis)
% torque 3 (along y-axis)

%% Constant Definition
m1 = 0.8; %kg
m3 = 68; %kg
l1 = 0.263; %m
l3 = 0.88; %m

b1 = 0.5;
b2 = 0.5;
b3 = 0.5;
%EPS = 1E-5;     % finite difference epsilon
% State Vector Derivative
xddot = AngleAccelN(b1,b2,b3,x(4,:),x(5,:),x(6,:),l1,l3,m1,m3,x(1,:),x(2,:),...
    x(3,:),u(1,:),u(2,:),u(3,:));
xdot = [x(4:6,:);xddot];

%Partial Derivatives
if nargout>1
%     x1 = repmat(x, [1,6]) + eye(6)*EPS;
%     x2 = repmat(x, [1,6]) - eye(6)*EPS;
%     uu = repmat(u,[1,6]);
%     f1 = leg_dynN(x1,uu);
%     f2 = leg_dynN(x2,uu);
%     xdot_x = (f1-f2)/2/EPS;
%     
%     xdot_u = zeros(6,3);
%     xdot_u(4:6,:) = eye(3);
      xdot_x = jacob_State(b1,b2,b3,x(4,:),x(5,:),x(6,:),l1,l3,...
          m1,m3,x(1,:),x(2,:),x(3,:),u(1,:),u(2,:),u(3,:));
      xdot_u = jacob_Cont(l1,l3,m1,m3,x(1,:),x(2,:),x(3,:));
end