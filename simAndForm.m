R1 = 150;
R2 = 0.4*R1;
lambda = 0.001;

theta = 15*pi/180; %0; %acos( 1 - R2*R2/(2*R1*R1) );
x = R1*cos(theta)-sqrt(R2*R2 - R1*R1*sin(theta)^2):3:R1;

probability_formula = zeros(1,numel(x));

for j = 1:numel(x)
    r = x(j) ;

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
    f1 = @(alpha) (1/2/pi)*exp(-pi*lambda*r11(alpha).*r11(alpha));
    f2 = @(alpha) (1/2/pi)*exp(-pi*lambda*r12(alpha).*r12(alpha));
    Prob = integral(f1,-a2,a1) +  integral(f2,a1,2*pi-a2);
    probability_formula(j) = Prob;
end
plot(probability_formula,'r')
hold on
%% Simulation

origin = [0,0];
L = 2000; 
area = L*L;

center1 = origin;
center2 = [R1,0];
%x = R1-R2:5:R1; % along the line joining the centers

probability_sim = zeros(1,numel(x));

for j = 1:numel(x)
     
    r = x(j);
    currentPosition = [r*cos(theta), r*sin(theta)];

    noOfIters = 4000;

    out = zeros(1,noOfIters);


    for i = 1:noOfIters

        N = poissrnd(lambda*area); % no. of AUs
        p = unifrnd(-L/2,L/2,N,2);

        distances = sqrt(sum((p-repmat(currentPosition,length(p),1))'.^2));
        minDist = min(distances);
        nextPosIndex = find(distances == minDist);
        nextPosition = p(nextPosIndex,:);

        if (norm(nextPosition-center1) > R1) || (norm(nextPosition-center2) > R2)
            out(i) = 1;
        end

    end

probability_sim(j) = sum(out)/noOfIters;

end
plot(probability_sim)