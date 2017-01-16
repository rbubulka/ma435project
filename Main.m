clc
clear all
close all

%%%%%%%%%%%%%%%%
% Functions

f = @(x) (1/2)*x^2;

MM = @(x,y) (1/2)*(sign(x)+sign(y))*min(abs(x),abs(y));

%vp = @(v1,v2) 0;
vp = @(v1,v2) MM(v1,v2);
%fp = @(x,y) 0;
fp = @(f1,f2) MM(f1,f2);

%%%%%%%%%%%%%%%%%
% Definitions
dt = .005;
dx = .005;

lambda = dt/dx;
tf = 0.8;

x1 = -2;
x2 = 2;

x = x1;
t = 0;

tstep = tf/dt;
xstep = (x2-x1)/dx;

v=zeros(xstep*2+1,tstep*2+2);
%%%%%%%%%%%%%%%%%
% Initialization

for j = 1:1:2*xstep
    if x<0
        v(j,1)=1;
    elseif x>0
        v(j,1)=0;
    end
%     v(j,1)=sin(pi*x);
    x=x+dx/2;
end

%%%%%%%%%%%%%%%%%

x=x1;
for i = 1:2:2*tstep
    for j = 3:2:2*xstep
        v(2*xstep+1,i)=0;
        %v(2*xstep,i);
        v(j,i+1)=v(j,i)-(1/2)*lambda*fp(v(j+2,i)-v(j,i),v(j,i)-v(j-2,i));
        v(j,i+2)=v(j,i)-(1/2)*lambda*fp(v(j+2,i)-v(j,i),v(j,i)-v(j-2,i));
        T1=(1/2)*(v(j,i)+v(j+2,i));
        
        T2=(1/8)*(vp(v(j+2,i)-v(j,i),v(j,i)-v(j-2,i))-vp(v(j+2,i)-v(j,i),v(j,i)-v(j-2,i)));
        
        
        T3=-lambda*(f(v(j+1,i+1))-f(v(j,i+1)));

        v(j+1,i+2)=T1+T2+T3;
        v(j+1,i+3)=T1+T2+T3;
        x=x+dx;
    end
    t=t+dt;
end

figure(1)
plot(linspace(0,2,j+2),v(:,11),'ob')

figure(2)
t=0
for i=1:size(v,2)-2
    plot(linspace(x1,x2,length(v)),v(:,i),'ok')
    text=sprintf('Time %2.3f',t);
    title(text)
    axis([x1,2,-4,4])
    pause(.25)
    t=t+dt/2;
end

