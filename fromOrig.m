L = 100; 
area = L*L;

lambda = 0.01; % rwp density

currentPosition = [0,0];
origin = [0,0];
r0 = 10;

distTravelled = 0;

while distance(currentPosition, origin) < r0
    
    N = poissrnd(lambda*area); % no. of AUs
    p = unifrnd(-L/2,L/2,N,2);

    distances = distance(p,currentPosition);
    minDist = min(distances);
    nextPosIndex = find(distances == minDist);
    nextPosition = p(nextPosIndex,:);
    
    currentPosition = nextPosition;
    distTravelled = distTravelled + minDist;
    
end

distTravelled

