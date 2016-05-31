function r0 = incircleLength(i1,i2,radius,center)
%% 
% given initial position i1 and final position i2
% i1 = [1,0];
% i2 = [-1,2];
    i1 = i1-center;
    i2 = i2-center;
    if norm(i2) < radius
        r0 = norm(i2-i1);
    else
        r2 = norm(i1);
        if r2 == 0
            r0 = radius;
        else
        r1 = radius; %radius
        theta=acos(dot(i2-i1,i1)/(r2*norm(i2-i1)));

        r0 = -r2*cos(theta) + sqrt( r1^2 - (r2*sin(theta))^2);
        end
    end
end
