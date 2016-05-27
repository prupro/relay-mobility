origin = [0,0];
L = 100; 
area = L*L;
lambda = 0.01; % rwp density

noOfLegs = 100;

noOfPaths = 10;
totalDist = zeros(1,noOfPaths);

for j = 1:noOfPaths

distTravelled = zeros(1, noOfLegs);

currentPosition = origin;

nextPosition = zeros( noOfLegs,2);

for i = 1:noOfLegs
    
    N = poissrnd(lambda*area); % no. of AUs
    p = unifrnd(-L/2,L/2,N,2);

    distances = sqrt(sum((p-repmat(currentPosition,length(p),1))'.^2));
    minDist = min(distances);
    nextPosIndex = find(distances == minDist);
    nextPosition(i,:) = p(nextPosIndex,:);

    currentPosition = nextPosition(i,:);
    distTravelled(i) = distTravelled(i) + minDist;
   
end

totalDist(j) = sum(distTravelled);
end
sum(totalDist)/(noOfLegs*noOfPaths)