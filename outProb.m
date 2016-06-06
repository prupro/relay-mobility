R1 = 10;
R2 = R1;
lambda = 0.01;

beta1 = @(theta) acos( 1 - R2*R2/(2*R1*R1) ) - theta ;
beta2 = acos( 1 - R2*R2/(2*R1*R1) );

alpha1 = @(r,theta) theta + atan( R1*sin(beta1(theta))/(R1*cos(beta1(theta)) - r) );
alpha2 = @(r,theta) atan( (r*sin(theta) + R1*sin(beta2))/(R1*cos(beta2) - r*cos(theta)));

r11 = @(r,theta,alpha) -r*cos(theta-alpha) + sqrt( R1*R1 - r*r * sin(theta-alpha)*sin(theta-alpha) );
r12 = @(r,theta,alpha) R1*cos(alpha) - r*cos(alpha-theta) + sqrt( (R1*cos(alpha)-r*cos(alpha-theta))*(R1*cos(alpha)-r*cos(alpha-theta)) - r*r*sin(theta)*sin(theta) - (R1-r*cos(theta))^2 + R2*R2 );

a1 = alpha1(r,theta);
a2 = alpha2(r,theta);
if a1 < 0
    a1 = pi + a1;
    a2 = pi + a2;
end
f1 = @(alpha) (1/2*pi)*exp(-2*pi*lambda*r11(r,theta,alpha)*r11(r,theta,alpha));
f2 = @(alpha) (1/2*pi)*exp(-2*pi*lambda*r12(r,theta,alpha)*r12(r,theta,alpha));
Prob = integral(@(alpha)f1,-a2,a1) +  integral(@(alpha)f2,2*pi-a2);
