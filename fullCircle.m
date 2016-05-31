origin = [0,0];
L = 1000; 
area = L*L;

lambda = 0.001; % rwp density

radius = 10;

center = origin;


noOfIter = 1000;
distTravelled = zeros(1, noOfIter);
noOfLegs = zeros(1,noOfIter);

for i = 1:noOfIter
    
    currentPosition = center;
    
    while (norm(currentPosition-center)<=radius)

        N = poissrnd(lambda*area); % no. of AUs
        p = unifrnd(-L/2,L/2,N,2);
        noOfLegs(i) = noOfLegs(i)+1;

        distances = sqrt(sum((p-repmat(currentPosition,length(p),1))'.^2));
        minDist = min(distances);
        nextPosIndex = find(distances == minDist);
        %currentPosition
        nextPosition = p(nextPosIndex,:);

        prevPosition = currentPosition;
        currentPosition = nextPosition;
        
        distTravelled(i) = distTravelled(i) + minDist;

    end
    distTravelled(i) = distTravelled(i) - minDist + incircleLength(prevPosition, currentPosition, radius, center);
end

averageDist = sum(distTravelled)/noOfIter
theoriticalDist = (1/(2*sqrt(lambda)))*(1- 2*qfunc(radius*sqrt(2*pi*lambda)))
