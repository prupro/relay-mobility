origin = [0,0];
L = 100; 
area = L*L;

lambda = 0.01; % rwp density

noOfIter = 1000;

distTravelled = zeros(1, noOfIter);

currentPosition = origin;

nextPosition = zeros( noOfIter,2);

for i = 1:noOfIter
    
    N = poissrnd(lambda*area); % no. of AUs
    p = unifrnd(-L/2,L/2,N,2);

    distances = distance(p,currentPosition);
    minDist = min(distances);
    nextPosIndex = find(distances == minDist);
    nextPosition(i,:) = p(nextPosIndex,:);

    currentPosition = nextPosition(i,:);
    distTravelled(i) = distTravelled(i) + minDist;
   
end

averageDist = cumsum(distTravelled)./(1:length(distTravelled));
plot(averageDist)