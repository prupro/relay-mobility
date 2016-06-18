origin = [0,0];
L = 1000; 
area = L*L;

R1 = 150;
R2 = 0.4*R1;
lambda = 0.001;


center1 = origin;
center2 = [R1,0];

theta = 5*pi/180;
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
    end
    if a2 < 0
        a2 = pi + a2;
    end
    f1 = @(alpha) (1/2/pi)*exp(-pi*lambda*r11(alpha).*r11(alpha));
    f2 = @(alpha) (1/2/pi)*exp(-pi*lambda*r12(alpha).*r12(alpha));
    Prob = integral(f1,-a2,a1) +  integral(f2,a1,2*pi-a2);
    probability_formula(j) = Prob;
end


noOfIter = 1000;

noOfLegs = zeros(1,numel(x));

for j = 1:numel(x)

    for i = 1:noOfIter

        currentPosition = [x(j)*cos(theta),x(j)*sin(theta)];

        while (norm(currentPosition-center1)<=R1)&&(norm(currentPosition-center2)<=R2+0.001)

            N = poissrnd(lambda*area); % no. of AUs
            p = unifrnd(-L/2,L/2,N,2);
            noOfLegs(j) = noOfLegs(j)+1;

            distances = sqrt(sum((p-repmat(currentPosition,length(p),1))'.^2));
            minDist = min(distances);
            nextPosIndex = find(distances == minDist);
            %currentPosition
            nextPosition = p(nextPosIndex,:);

            prevPosition = currentPosition;
            currentPosition = nextPosition;      

        end

    end
end
noOfLegs = noOfLegs/noOfIter;

plot(x,1./probability_formula,'r',x,noOfLegs);
legend('analytic','simulation');
title('No. of steps to leave;\theta=5^0, \lambda=0.001');
xlabel('Distance from BS');
ylabel('# steps');

% save('figs/oneOutProbl0.0005t15.txt', 'x', '-ASCII','-append');
% save('figs/oneOutProbl0.0005t15.txt', 'probability_formula', '-ASCII','-append');
% save('figs/oneOutProbl0.0005t15.txt', 'probability_sim', '-ASCII','-append');


