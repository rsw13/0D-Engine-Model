function lift = cam_profile(theta)

theta = deg2rad(theta);
r1 = 2;
r2 = 2;
cc = 1;

a = asin((r1 - r2)/cc);

if theta > pi/2
    
    beta = asin(r1 / (cc + r2));
    lift = (cc + r2) * cos(beta + (pi / 2) - theta) - r1;
    
    if lift < 0

        lift = 0;        

    end
    
    return
end

lamda = (tan(theta) * r1 * cos(a) - r1 * sin(a)) / (1/(cos(a)) +...
    tan(theta) * tan(a));

s = [sin(a); cos(a)] + lamda * [(1/cos(a)); -(sin(a)) / (cos(a))];

    if s(1) > cc

       l = cc ;
       y = cc / tan(theta);

       c = [1, (2 * y * cos(theta)), (y^2 - r2^2)]; 

       d(1) = (-c(2) + sqrt(c(2)^2 - 4 * c(1) * c(3))) / ( 2 * c(1));
       d(2) = (-c(2) - sqrt(c(2)^2 - 4 * c(1) * c(3))) / ( 2 * c(1));

       if isreal(d(1)) && (d(1) >= 0)

           d = d(1);

       elseif isreal(d(2)) && (d(2) >= 0)

           d= d(2);

       end

       s(1) = l + d * sin(theta);
       s(2) = y + d * cos(theta);

    end

    lift = sqrt( (s(1)^2 + s(2) ^2) - r1 ^ 2) - r1;

    if lift < 0

        lift = 0;        

    end

end

