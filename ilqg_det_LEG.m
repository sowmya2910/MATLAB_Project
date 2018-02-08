%[x_, u_, L, cost ] = 
%       ilqg_det(fnDyn, fnCost, dt, n, x0, u0, uMin, uMax)
%
%Iterative LQG design for nonlinear plants with bounded controls
%
%DYNAMICS:  x(t+1) = x(t) + dt*fnDyn(x(t), u(t))
%               where x(1) = x0, t = 1:n-1
%COST:      fnCost(x(n)) + sum_t dt*fnCost(x(t),u(t),t)
%CONTROL:   u(t) = u_(t) + L(:,:,t) * (x(t) - x_(t))
%
%INPUTS:
% fnDyn - handle of dynamics function, which must be in the format:
%           [xdot, xdot_x, xdot_u] = fnDyn(x, u)
%          where x is state, u is control, xdot is time derivative of x,
%          xdot_x and xdot_u are partial derivatives of xdot
%          if only xdot is required, x and u may be matrices
% fnCost- handle of cost function, which must be in the format:
%           [l, l_x, l_xx, l_u, l_uu, l_ux] = fnCost(x, u, t)
%          where l is the cost, and l_x etc. are its partial derivatives
%          t is the current time step; t=NaN means final cost
%          if only l is required, x and u may be matrices
% dt    - time step for discrete-time simulation
% n     - number of time steps
% x0    - initial state
% u0    - initial control, either vector or matrix with n-1 columns; if u0
%          is vector, it is repeated n-1 times to form a control sequence
% uMin  - [optional] minimum control, scalar or vector (-Inf: unbounded)
%          if scalar, the bound is the same for all elements of u
% uMax  - [optional] maximum control, scalar or vector (+Inf: unbounded)
%
%OUTPUTS:
% x_    - average trajectory under optimal controller (n column matrix)
% u_    - open-loop component of optimal controller (n-1 column matrix)
% L     - time-varying feedback gains of optimal controler (3D array)
% cost  - expected cost under optimal controller
%
%NOTES:
% A number of user-adjustable parameters controlling the iterative
%  algorithm are defined at the beginning of the file
%
% This is an implementation of the algorithm described in:
%  Todorov, E. and Li, W. (2005) A generalized iterative LQG method for
%  locally-optimal feedback control of constrained nonlinear stochastic
%  systems. In proceedings of the American Control Conference, pp 300-306
% The paper is available online at www.cogsci.ucsd.edu/~todorov
%
% THIS CODE ASSUMES A DETERMINISTIC SYSTEM

% Copyright (C) Emanuel Todorov, 2005-2006


function [x, u, L, cost] = ilqg_det_LEG(fnDyn, fnCost, dt, n, x0, u0, varargin)
global flgFix;

%---------------------- user-adjustable parameters ------------------------

lambdaInit = 0.01;       % initial value of Levenberg-Marquardt lambda
lambdaFactor = 5;      % factor for multiplying or dividing lambda
lambdaMax = 10000;       % exit if lambda exceeds threshold
epsConverge = 0.001;    % exit if relative improvement below threshold
maxIter = 1000;          % exit if number of iterations reaches threshold
flgPrint = 2;           % show cost- 0:never, 1:every iter, 2:final
maxValue = 1E80;       % upper bound on states and costs (Inf: none)


%---------------------------------------------- get optional arguments
if nargin>=7,
    uMin = varargin{1};
else
    uMin = -Inf;
end
if nargin>=8,
    uMax = varargin{2};
else
    uMax = Inf;
end


%---------------------------------------------- initialization
szX = size(x0, 1);          % size of state vector
szU = size(u0, 1);          % size of control vector

L = zeros(szU, szX, uint64(n-1));   % init feedback gains

if size(u0,2)==1,           % init control sequence
    u = repmat(u0, [1,n-1]);
else
    u = u0;
end

if isscalar(uMin),          % transform scalar arguments into vectors
    uMin = repmat(uMin, [szU,1]);
end
if isscalar(uMax),
    uMax = repmat(uMax, [szU,1]);
end

flgFix = 0;                 % clear large-fix global flag

x = zeros(szX, n);          % init state sequence and cost
[x, cost] = simulate(fnDyn, fnCost, dt, x0, u, maxValue);


%---------------------------------------------- optimization loop
lambda = lambdaInit;
flgChange = 1;

for iter = 1:maxIter,
    
    %------ STEP 1: approximate dynamics and cost along new trajectory
    if flgChange,
        [s0(n),s(:,n),S(:,:,n)] = fnCost(x(:,n), NaN, NaN);  % final cost
        
        for k = n-1:-1:1
            % quadratize cost, adjust for dt
            [l0,l_x,l_xx,l_u,l_uu,l_ux] = fnCost(x(:,k), u(:,k), k);
            q0(k) = dt * l0;
            q(:,k) = dt * l_x;
            Q(:,:,k) = dt * l_xx;
            r(:,k) = dt * l_u;
            R(:,:,k) = dt * l_uu;
            P(:,:,k) = dt * l_ux;

            % linearize dynamics, adjust for dt
            [f, f_x, f_u] = fnDyn(x(:,k), u(:,k));
            A(:,:,k) = eye(szX) + dt * f_x;
            B(:,:,k) = dt * f_u;
        end

        flgChange = 0;
    end

    
    %------ STEP 2: compute optimal control law and cost-to-go
    for k = n-1:-1:1
        % compute shortcuts g,G,H
        g = r(:,k) + B(:,:,k)'*s(:,k+1);
        G = P(:,:,k) + B(:,:,k)'*S(:,:,k+1)*A(:,:,k);
        H = R(:,:,k) + B(:,:,k)'*S(:,:,k+1)*B(:,:,k);

        % find control law
        [l(:,k), L(:,:,k)] = uOptimal(g,G,H,u(:,k),uMin,uMax,lambda);

        % update cost-to-go approximation
        S(:,:,k) = Q(:,:,k) + A(:,:,k)'*S(:,:,k+1)*A(:,:,k) + ...
            L(:,:,k)'*H*L(:,:,k) + L(:,:,k)'*G + G'*L(:,:,k);
        s(:,k) = q(:,k) + A(:,:,k)'*s(:,k+1) + ...
            L(:,:,k)'*H*l(:,k) + L(:,:,k)'*g + G'*l(:,k);
        s0(k) = q0(k) + s0(k+1) + l(:,k)'*H*l(:,k)/2 + l(:,k)'*g;

        % HACK USED TO PREVENT OCCASIONAL DIVERGENCE
        if isfinite(maxValue),
            S(:,:,k) = fixBig(S(:,:,k), maxValue);
            s(:,k) = fixBig(s(:,k), maxValue);
            s0(k) = fixBig(s0(k), maxValue);
        end
    end

    
    %------ STEP 3: new control sequence, trajectory, cost
    % simulate linearized system to compute new control
    dx = zeros(szX,1);
    for k=1:n-1
        du = l(:,k) + L(:,:,k)*dx;   
        du = min(max(du+u(:,k),uMin),uMax) - u(:,k);
        dx = A(:,:,k)*dx + B(:,:,k)*du;  
        unew(:,k) = u(:,k) + du;
    end

    % simulate system to compute new trajectory and cost
    [xnew, costnew] = simulate(fnDyn, fnCost, dt, x0, unew, maxValue);
    
    
    %------ STEP 4: Levenberg-Marquardt method
    if costnew<cost,
        % decrease lambda (get closer to Newton method)
        lambda = lambda / lambdaFactor;
        
        % accept changes, flag changes
        u = unew;
        x = xnew;
        flgChange = 1;
        
        if flgPrint==1,
            fprintf('Iteration = %d;  Cost = %.4f;  logLambda = %.1f\n', ...
                iter, costnew, log10(lambda) );
        end

        if iter>1 & (abs(costnew - cost)/cost < epsConverge),
            cost = costnew;
            break;            % improvement too small: EXIT
        end
        cost = costnew;

    else
        % increase lambda (get closer to gradient descent)
        lambda = lambda * lambdaFactor;
        
        if lambda>lambdaMax
            break;            % lambda too large: EXIT
        end
    end
end


% print final result if necessary
if flgPrint==2,
    fprintf('Iterations = %d;  Cost = %.4f\n', iter, cost);
end

if flgFix,
    warning('ilqg_det had to limit large numbers, results may be inaccurate');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Simulate controlled system, compute trajectory and cost
function [x, cost] = simulate(fnDyn, fnCost, dt, x0, u, maxValue)
endSim = 0;
minHeight = 0.25*sqrt(2)/2;
% get sizes
szX = size(x0,1);
szU = size(u,1);
n = size(u,2)+1;

% initialize simulation
x = zeros(szX, n);
x(:,1) = x0;
cost = 0;
% Bound determination
Ry = @(q1) [[cos(q1) 0 sin(q1)];[0 1 0];[-sin(q1) 0 cos(q1)]];
Rx = @(q3) [[1 0 0];[0 cos(q3) -sin(q3)];[0 sin(q3) cos(q3)]];
Rleg = @(q1,q2,q3) Ry(q1) * Rx(q2) * Ry(q3);
l1 = 0.263;
l3 = 0.88;
e1 = l1/2 *[-1;0;0];
e2 = l1 * [-1;0;0];
e3 = l3*[0;0;-1];

% run simulation with substeps
for k = 1:n-1
    if ~endSim
        x(:,k+1) = x(:,k) + dt * fnDyn(x(:,k), u(:,k));
        revs = x(1:3,k+1)./(2*pi);
        frac = mod(revs,sign(revs));
        for n = 1:3
            if frac(n)~=0
                x(n,k+1) = frac(n).*2*pi;
            end
        end

    else
        x(:,k+1) = x(:,k);
        revs = x(1:3,k+1)./(2*pi);
        frac = mod(revs,1);
        x(1:3,k+1) = frac.*2*pi;
    end
    q1 = x(1,k+1);
    q2 = x(2,k+1);
    q3 = x(3,k+1);
    p2 = Ry(q1)*e2;
    p3 = p2+Rleg(q1,q2,q3)*e3;
    if p3(3) > minHeight || q1>pi/1.7 || q1<-pi/2
        endSim = 1;
    end
    x(:,k+1) = fixBig(x(:,k+1), maxValue);
    cost = cost + dt * fnCost(x(:,k), u(:,k), k);
end

% adjust for final cost, subsample trajectory
cost = cost + fnCost(x(:,end), NaN, NaN);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Compute optimal control law
function [l,L] = uOptimal(g, G, H, u, uMin, uMax, lambda)

% eigenvalue decomposition, modify eigenvalues
[V,D] = eig(H);
d = diag(D);
d(d<0) = 0;
d = d + lambda;

% inverse modified Hessian, unconstrained control law
H1 = V*diag(1./d)*V';
l = -H1*g;
L = -H1*G;

% enforce constraints
l = min(max(l+u,uMin),uMax) - u;

% modify L to reflect active constraints
L((l+u==uMin)|(l+u==uMax),:) = 0;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Limit numbers that have become suspiciously large
function s = fixBig(s, maxValue)
global flgFix;

ibad = (abs(s)>maxValue);
s(ibad) = sign(s(ibad)) * maxValue;

if any(ibad(:)),
    flgFix = 1;
end
