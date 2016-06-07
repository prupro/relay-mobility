origin = [0,0];
L = 2000; 
area = L*L;
lambda = 0.001;

R1 = 150;
R2 = 0.4*R1;
center1 = origin;
center2 = [R1,0];
%x = R1-R2:5:R1; % along the line joining the centers

theta = 15*pi/180;
x = R1*cos(theta)-sqrt(R2*R2 - R1*R1*(sin(theta))^2):3:R1;
%r = R1 - R2/2;
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


