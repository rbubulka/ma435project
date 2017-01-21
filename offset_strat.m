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
dt = .002;
dx = .001;

lambda = dt/dx;
tf = 0.8;

x1 = -2;
x2 = 2;

offset = 0;
x = x1-2*dx;
t = 0;

tstep = tf/dt;
xstep = (x2-x1)/dx;

v=zeros(xstep*2+1,tstep*2+4);
%%%%%%%%%%%%%%%%%
% Initialization

for j = 1:1:2*xstep+4
    if x<0
        v(j,1)=1;
    elseif x>0
        v(j,1)=0;
    end
%     v(j,1)=sin(pi*x);
    x=x+dx/2;
end
v(1,2:1:size(v,2))=1;
v(2,2:1:size(v,2))=1;


% v(:,2)=v(:,1);

%%%%%%%%%%%%%%%%%

clear x
for i = 1:2:2*tstep
    for j = 3:2:2*xstep
        v(j+offset,i+1)=v(j+offset,i)-(1/2)*(lambda)*fp(f(v(j+2+offset,i))-f(v(j+offset,i)),f(v(j+offset,i))-f(v(j-2+offset,i)));
        
        T1=(1/2)*(v(j+offset,i)+v(j+2+offset,i));
        T2=(1/8)*(vp(v(j+2+offset,i)-v(j+offset,i),v(j+offset,i)-v(j-2+offset,i))-vp(v(j+4+offset,i)-v(j+2+offset,i),v(j+2+offset,i)-v(j+offset,i)));
        T3=-lambda*(f(v(j+2+offset,i+1))-f(v(j+offset,i+1)));
        
        TT1(j,i)=T1;
        TT2(j,i)=T2;
        TT3(j,i)=T3;

        v(j+1+offset,i+2)=T1+T2+T3;
    end
    v(2,i+2)=v(4,i+2);
    v(2,i+3)=v(4,i+2);
    v(3,i)=v(5,i+1);
    v(3,i+1)=v(5,i+1);
    if offset == 0
        offset = 1;
    elseif offset == 1
        offset = 0;
    end
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