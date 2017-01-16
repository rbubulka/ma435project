% this is a test

% test commit

%%%%%%%%%%%%%%%%
% Functions

f = @(x) (1/2)*x^2;

MM = @(x,y) (1/2)*(sign(x)+sign(y))*min(x,y);

%vp = @(v(x+1/2,t),v(x-1/2,t))
vp = @(v1,v2) MM(v1,v2);
fp = @(f1,f2) MM(f1,f2);

%%%%%%%%%%%%%%%%%
% Definitions
dt = .001;
dx = .001;

lambda = dt/dx;
tf = 0.8;

x = 0;
t = 0;

x1 = 0;
x2 = 2;

tstep = tf/dt;
xstep = (x2-x1)/dx;
%%%%%%%%%%%%%%%%%
% Initialization

for j = 1:2:(xstep+4)/2
    v(j,1)=sin(pi*x);
    x=x+dx/2;
end

%%%%%%%%%%%%%%%%%

x=0;
for i = 1:2:tstep/2
    for j = 3:2:xstep/2
        v(j,i+1)=v(j,i)-(1/2)*lambda*fp(j,i);
        T1=(1/2)*(v(j,i)+v(j+2,i));
        T2=(1/8)*(vp(x,t)-vp(x+dt,t));
        T3=-lambda*(f(v(j+1,i+1))-f(v(j,i+1)));

        v(j+1,i+2)=T1+T2+T3;
        x=x+dx;
    end
    t=t+dt;
end

