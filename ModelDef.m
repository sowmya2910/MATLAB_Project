%% Rotation Matricies
g = 9.81;
%m1 = 1.2;
%m3 = 68;
%ft_l = 0.2;
%leg_l = 0.88;
syms l1 l3 m1 m3

assume(l1,'real')
assumeAlso(l1>0)
assume(l3,'real')
assumeAlso(l3>0)
assume(m1,'real')
assumeAlso(m1>0)
assume(m3,'real')
assumeAlso(m3>0)




Ry = @(q1) [[cos(q1) 0 sin(q1)];[0 1 0];[-sin(q1) 0 cos(q1)]];
Rx = @(q3) [[1 0 0];[0 cos(q3) -sin(q3)];[0 sin(q3) cos(q3)]];
Rleg = @(q1,q2,q3) Ry(q1) * Rx(q2) * Ry(q3);
S = @(x,y,z) [[0 -z y];[z 0 -x];[-y x 0]];


%% Ex
syms th1 th2 th3 dth1 dth2 dth3 ddth1 ddth2 ddth3

assume(th1,'real')
assume(th2,'real')
assume(th3,'real')
assume(dth1,'real')
assume(dth2,'real')
assume(dth3,'real')
assume(ddth1,'real')
assume(ddth2,'real')
assume(ddth3,'real')

e1 = l1/2 *[-1;0;0];
e2 = l1 * [-1;0;0];
e3 = l3*[0;0;-1];
p1 = Ry(th1)*e1;
p2 = Ry(th1)*e2;
p3 = p2+Rleg(th1,th2,th3)*e3;


dp1 = S(0,dth1,0)*Ry(th1)*e1;
dp2 = S(0,dth1,0)*Ry(th1)*e2;
J = @(th1,th2,th3) [[sin(th1)*sin(th2) cos(th3) 0];[cos(th2) 0 1];...
    [-cos(th1)*sin(th2) sin(th3) 0]];
om = J(th1,th2,th3) * [dth1;dth2;dth3];
dp3 = dp2+Rleg(th1,th2,th3)*cross(om,e3);

T = 0.5*dp1'*m1*dp1 + 0.5*dp3'*m3*dp3;
V = [0 0 -1]*m1*g*p1 + ...
    [0 0 -1]*m3*g*p3;
L = T-V;

ddq = [ddth1;ddth2;ddth3];
dq = [dth1;dth2;dth3];
q = [th1;th2;th3];
tor = [th1;th2;th3];
for i = 1:3
    dqiL = diff(L,dq(i));
    ddtL(i,1) = dth1*diff(dqiL,th1)+dth2*diff(dqiL,th2)+dth3*diff(dqiL,th3) + ...
        ddth1*diff(dqiL,dth1)+ddth2*diff(dqiL,dth2)+ddth3*diff(dqiL,dth3);
    dqL(i,1) = diff(L,q(i));
    tor(i,1) = ddtL(i,1) - diff(L,q(i));
end
    
M = simplify(jacobian(ddtL,[ddth1,ddth2,ddth3]));

for i = 1:3
    for j = 1:3
        ent = 0;
        for k = 1:3
            ent = ent+(diff(M(i,j),q(k))+diff(M(i,k),q(j))-diff(M(k,j),q(i)))*dq(k);
        end
        ent = 0.5*ent;
        if ~isempty(ent)
            C(i,j) = ent;
        else
            C(i,j) = 0;
        end
    end
end

for i = 1:3
    G(i,1) = diff(V,q(i));
end
syms b1 b2 b3
assume(b1,'real')
assume(b2,'real')
assume(b3,'real')

N = [b1*dq(1);b2*dq(2);b3*dq(3)];
%u_mod = M*ddq+C*dq+N+G;
syms u1 u2 u3
assume(u1,'real')
assume(u2,'real')
assume(u3,'real')

ddq_mod = simplify(inv(M)*([u1;u2;u3]-C*dq-G));
ddq_modN = simplify(inv(M)*([u1;u2;u3]-C*dq-N-G));
% matlabFunction(ddq_mod,'File','AngleAccel')
% matlabFunction(ddq_modN,'File','AngleAccelN')

% after rebuilding AngleAccel functions,
% make sure to change input order too

%% Partial Derivatives
bigF = [dth1;dth2;dth3;ddq_modN];
dbigF_x = simplify(jacobian(bigF,[th1;th2;th3;dth1;dth2;dth3]));
dbigF_u = simplify(jacobian(bigF,[u1;u2;u3]));
matlabFunction(dbigF_x,'File','jacob_State')
matlabFunction(dbigF_u,'File','jacob_Cont')