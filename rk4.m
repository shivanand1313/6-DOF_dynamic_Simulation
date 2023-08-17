
function [t,y] = rk4(dydt, tf, dt, y0)
% reshape y0 to column vector
y0 = reshape(y0, [], 1);

n = tf/dt;
d = size(y0,1);

t = zeros(n,1);
y = zeros(n,d);

tn = 0;
yn = y0;

t(1,:) = tn;
y(1,:) = yn';

for i=2:n+1
C1 = dydt(tn, yn);
C2 = dydt((tn + dt*0.5), (yn + dt*C1*0.5));
C3 = dydt((tn + dt*0.5), yn + (dt*C2*0.5));
C4 = dydt(tn + dt, yn + dt*C3);

yn = yn + (C1 + 2*C2 + 2*C3 + C4)*dt/6;
tn = tn + dt;

t(i) = tn;
y(i,:) = yn';
end
end
