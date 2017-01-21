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
dx = .01;

lambda = dt/dx;
tf = 0.8;

x1 = -.1;
x2 = .1;

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
v(:,2)=v(:,1);

%%%%%%%%%%%%%%%%%

clear x
for i = 1:1:2*tstep
    for j = 3:1:2*xstep-3
        v(j,i+1)=v(j,i)-(1/2)*lambda*fp(f(v(j+2,i))-f(v(j,i)),f(v(j,i))-f(v(j-2,i)));
%         v(j,i+2)=v(j,i)-(1/2)*lambda*fp(f(v(j+2,i))-f(v(j,i)),f(v(j,i))-f(v(j-2,i)));
        
        T1=(1/2)*(v(j,i)+v(j+2,i));
        T2=(1/8)*(vp(v(j+2,i)-v(j,i),v(j,i)-v(j-2,i))-vp(v(j+4,i)-v(j+2,i),v(j+2,i)-v(j,i)));
        T3=-lambda*(f(v(j+2,i+1))-f(v(j,i+1)));

        v(j+1,i+2)=T1+T2+T3;
%         v(j+1,i+3)=T1+T2+T3;
    end
    v(2,i+2)=v(4,i+2);
    v(1,i+1)=v(3,i+1);
%     v(j+2,i+1)=v(j,i+1);
%     v(j+3,i+2)=v(j+1,i+2);
t=t+dt;
end

% figure(1)
% plot(linspace(0,2,j+4),v(:,11),'ob')

figure(2)
t=0
x=linspace(x1,x2,size(v,1));
for i=1:size(v,2)+2
    plot(x(1:2:length(x)),v(1:2:size(v,1),i),'ob',x(2:2:length(x)),v(2:2:size(v,1),i),'xr')
    text=sprintf('Time %2.3f',t);
    title(text)
    axis([x1,x2,-4,4])
    pause(.25)
    t=t+dt/2;
end