R1 = 150;
R2 = 0.4*R1;
r = (R1-R2)/2;
theta = 0; %acos( 1 - R2*R2/(2*R1*R1) );
lambda = 0.01;


beta1 = acos( 1 - R2*R2/(2*R1*R1) ) - theta;
beta2 = acos( 1 - R2*R2/(2*R1*R1) );


alpha1 = theta + atan( R1*sin(beta1)/(R1*cos(beta1) - r) );
alpha2 = atan( (r*sin(theta) + R1*sin(beta2))/(R1*cos(beta2) - r*cos(theta)));

r11 = @(alpha) -r*cos(theta-alpha) + sqrt( R1*R1 - r*r * sin(theta-alpha).*sin(theta-alpha) );
r12 = @(alpha) R1*cos(alpha) - r*cos(alpha-theta) + sqrt( (R1*cos(alpha)-r*cos(alpha-theta)).*(R1*cos(alpha)-r*cos(alpha-theta)) - r*r*sin(theta)*sin(theta) - (R1-r*cos(theta))^2 + R2*R2 );

a1 = alpha1;
a2 = alpha2;
if a1 < 0
    a1 = pi + a1;
    a2 = pi + a2;
end
f1 = @(alpha) (1/2/pi)*exp(-2*pi*lambda*r11(alpha).*r11(alpha));
f2 = @(alpha) (1/2/pi)*exp(-2*pi*lambda*r12(alpha).*r12(alpha));
Prob = integral(f1,-a2,a1) +  integral(f2,a1,2*pi-a2)
