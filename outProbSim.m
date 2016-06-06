origin = [0,0];
L = 2000; 
area = L*L;
lambda = 0.001;

R1 = 150;
R2 = 0.4*R1;
center1 = origin;
center2 = [R1,0];

r = R1 - R2/2 +25;
theta = 0; 
currentPosition = [r*cos(theta), r*sin(theta)];

noOfIters = 10000;

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
simOutProb = sum(out)/noOfIters


