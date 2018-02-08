%Cost function for bending leg with leg_dyn
%
% [l, l_x, l_xx, l_u, l_uu, l_ux] = leg_cost(x, u, t, target)
%
% t = nan signals final time

function [l, l_x, l_xx, l_u, l_uu, l_ux] = leg_cost(x, u, t, target)


wp = 1E+9;      % terminal position cost weight
wv = 1E+3;      % terminal velocity cost weight


load('costfilter.mat')

st_costvals = @(x) polyval(y1f,x(1))+polyval(y2f,x(2))+polyval(y3f,x(3));
% compute cost
if isnan(t)                        % final cost
    l = wp*sum((x(1:3)-target).^2) + wv*sum(x(4:6).^2);
else
    l = st_costvals(x)+sum(u.^2);
end


% compute derivatives of cost
if nargout>1

    l_x = zeros(6,1);
    l_x(1) = polyval(polyder(y1f),x(1));
    l_x(2) = polyval(polyder(y2f),x(2));
    l_x(3) = polyval(polyder(y3f),x(3));
    l_xx = zeros(6,6);
    l_xx(1,1) = polyval(polyder(polyder(y1f)),x(1));
    l_xx(2,2) = polyval(polyder(polyder(y2f)),x(2));
    l_xx(3,3) = polyval(polyder(polyder(y3f)),x(3));

    l_u = 2*u;
    l_uu = 2*eye(3);
    l_ux = zeros(3,6);
    
    if isnan(t)                           % final cost
        l_x(1:3) = 2*wp*(x(1:3)-target);
        l_x(4:6) = 2*wv*x(4:6);
        l_xx = 2*diag([wp wp wp wv wv wv]);
    end
end