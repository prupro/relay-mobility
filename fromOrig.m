origin = [0,0];
L = 100; 
area = L*L;

lambda = 0.01; % rwp density

noOfIter = 1000;

distTravelled = zeros(1, noOfIter);

for i = 1:noOfIter
    
    currentPosition = [0,0];

    r0 = 10;

    
    while distance(currentPosition, origin) < r0

        N = poissrnd(lambda*area); % no. of AUs
        p = unifrnd(-L/2,L/2,N,2);

        distances = distance(p,currentPosition);
        minDist = min(distances);
        nextPosIndex = find(distances == minDist);
        nextPosition = p(nextPosIndex,:);

        currentPosition = nextPosition;
        distTravelled(i) = distTravelled(i) + minDist;

    end
    
end

averageDist = cumsum(distTravelled)./(1:length(distTravelled));
plot(averageDist)
