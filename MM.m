function [w] = MM(x,y,z)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[j,k] = size(x);
for i=1:j
    for m=1:k
        if x(i,m) > 0 & y(i,m) > 0 & z(i,m) > 0
            input = [x(i,m),y(i,m),z(i,m)];
            w(i,m) = min(input);
        elseif x(i,m) < 0 & y(i,m) < 0 & z(i,m) < 0
            input = [x(i,m),y(i,m),z(i,m)];
            w(i,m) = max(input);
        else
            w(i,m) = 0;
        end

    end
end

end

