clc
clear all
close all

%%%%%%%%%%%%%%%%
% Functions

f = @(x) (1/2)*x.^2;
fp = @(x) (x);
g = @(y) (1/2)*y.^2;
gp = @(y) (y);

%w'jk = MM{delw(j-.5,k)+(1/2)*MM(del^2(w(j-1,k)),del^2(w(j,k))),delw(j+.5,k)-(1/2)*MM(del^2(w(j,k)),del^2(w(j+1,k)))}
%also need w'for y-.5.....


wpx = @(w,j,k,i) MM(w(j+1,k,i)-w(j,k,i),(1/2)*(w(j+1,k,i)-w(j-1,k,i)),w(j,k,i)-w(j-1,k,i));
wpy = @(w,j,k,i) MM(w(j,k+1,i)-w(j,k,i),(1/2)*(w(j,k+1,i)-w(j,k-1,i)),w(j,k,i)-w(j,k-1,i));
%%%%%%%%%%%%%%%%%
% Definitions
dt = .05;
dx = .1;
dy = .1;

lambda = dt/dx;
mew = dt/dy;
tf = 2;

x1 = -1;
x2 = 1;
y1 = -1;
y2 = 1;
z1 = -1;
z2 = 1;

offset = 0;
x = x1-2*dx;
t = 0;

tstep = tf/dt;
xstep = (x2-x1)/dx;
ystep = (y2-y1)/dy;

v=zeros(xstep*2+1,ystep*2+1,tstep*2+4);
%%%%%%%%%%%%%%%%%
% Initialization

for j = 1:1:2*xstep+4
    y=y1-2*dy;
    for k=1:1:2*xstep+4
        if x>0 & y>0
            v(j,k,1)=-1;
        elseif x<0 & y>0
            v(j,k,1)=-.2;
        elseif x<0 & y<0
            v(j,k,1)=.5;
        elseif x>0 & y<0
            v(j,k,1)=.8;
        end
        y=y+dy/2;
    end
    x=x+dx/2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BC

for i=3:2:size(v,2)
    v(1:2,:,i)=v(1:2,:,1);
    v(:,1:2,i)=v(:,1:2,1);
    v(43:44,:,i)=v(43:44,:,1);
    v(:,43:44,i)=v(:,43:44,1);
end

%%%%%%%%%%%%%%%%%

clear x
for i = 1:2:2*tstep
    for k = 3:2:2*ystep
        for j = 3:2:2*xstep
            
            I1 = v(j+offset,k+offset,i);
            I2 = -(1/2)*(lambda)*fp(v(j+offset,k+offset,i));
            I3 = -(1/2)*(mew)*gp(v(j+offset,k+offset,i));
            v(j+offset,k+offset,i+1) = I1+I2+I3;
            %fp(f(v(j+2+offset,i)-v(j+offset,i)),f(v(j+offset,i)-v(j-2+offset,i)));
        end
    end
    for k = 3:2:2*ystep
        for j = 3:2:2*xstep
            T1=(1/4)*(v(j+offset,k+offset,i)+v(j+1+offset,k+offset,i)+v(j+offset,k+1+offset,i)+v(j+1+offset,k+1+offset,i));
            T2=-(lambda/2)*(f(v(j+1+offset,k+offset,i+1))-f(v(j+offset,k+offset,i+1)));
            T3=-(lambda/2)*(f(v(j+1+offset,k+1+offset,i+1))-f(v(j+offset,k+1+offset,i+1)));
            T4=-(mew/2)*(g(v(j+offset,k+1+offset,i+1))-g(v(j+offset,k+offset,i+1)));
            T5=-(mew/2)*(g(v(j+1+offset,k+1+offset,i+1))-g(v(j+1+offset,k+offset,i+1)));
            T6=(1/16)*(wpx(v,j+offset,k+offset,i)-wpx(v,j+1+offset,k+offset,i)+wpx(v,j+offset,k+1+offset,i)-wpx(v,j+1+offset,k+1+offset,i));
            T7=(1/16)*(wpy(v,j+offset,k+offset,i)-wpy(v,j+offset,k+1+offset,i)+wpy(v,j+1+offset,k+offset,i)-wpy(v,j+1+offset,k+1+offset,i));

            v(j+1+offset,k+1+offset,i+2)=T1+T2+T3+T4+T5+T6+T7;
        end
    end
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
for i=1:size(v,3)+2
    surf(v(:,:,i))
    text=sprintf('Time %2.3f',t);
    title(text)
    axis([0,80,0,80,z1,z2])
    pause(.1)
    t=t+dt/2;
end