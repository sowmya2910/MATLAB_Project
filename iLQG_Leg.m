clear

initial_pos = [-pi/12;0;0];
initial_v = [-0.2;0;0];
target = [0;0;0] ;          % target configuration in radians

dt = 0.0001;                  % time step
t_dur = 0.2;
s_dur = 0.2;

bound_fitting
x0 = [initial_pos; initial_v];
u0 = [1;0;-1];
xStates = [];
uStates = [];
for k = 1:uint64(t_dur/s_dur)
    n = uint64(s_dur/dt);                     % number of time steps
    % form initial state and control

    fnCost = @(x_,u_,t_) leg_cost(x_,u_,t_,target);
    % solve optimization problem
    [x_out, u_out] = ilqg_det_LEG(@leg_dynN, fnCost, dt, n, x0, u0,-Inf,[0;Inf;Inf]);
    xStates = horzcat(xStates,x_out);
    uStates = horzcat(uStates,u_out,[0;0;0]);
    x0 = x_out(:,end);
end

%% Plots
figure(3)
plot(dt:dt:t_dur,uStates')
figure(4)

subplot(211)
plot(dt:dt:t_dur,xStates(1:3,:)')
subplot(212)
plot(dt:dt:t_dur,xStates(4:6,:)')

figure(5)
subplot(211)
plot(dt:dt:t_dur,xStates(1:3,:).*(180/pi)')
subplot(212)
plot(dt:dt:t_dur,xStates(4:6,:)')

%% XYZ Plots
Ry = @(q1) [[cos(q1) 0 sin(q1)];[0 1 0];[-sin(q1) 0 cos(q1)]];
Rx = @(q3) [[1 0 0];[0 cos(q3) -sin(q3)];[0 sin(q3) cos(q3)]];
Rleg = @(q1,q2,q3) Ry(q1) * Rx(q2) * Ry(q3);

l1 = 0.263;
l3 = 0.88;

e1 = l1/2 *[-1;0;0];
e2 = l1 * [-1;0;0];
e3 = l3*[0;0;-1];
% p1 = Ry(th1)*e1;
% p2 = Ry(th1)*e2;
% p3 = p2+Rleg(th1,th2,th3)*e3;
for n = 1:length(x_out(1,:))
    pos_ank(:,n) = Ry(x_out(1,n))*e2;
    pos_hip(:,n) = pos_ank(:,n)+Rleg(x_out(1,n),x_out(2,n),x_out(3,n))*e3;
end
figure(6)
scatter(0,0,15,'k*')
hold on
for i = 1:length(x_out(1,:))
%     h1 = line;
%     h2 = line;
%     h1.XData = [0 -pos_ank(1,i)];
%     h1.YData = [0 -pos_ank(3,i)];
%     h2.XData = [-pos_ank(1,i) -pos_hip(1,i)];
%     h2.YData = [-pos_ank(3,i) -pos_hip(3,i)];
%     pause(0.01)
      
end
hold off
    
    




