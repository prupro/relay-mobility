origin = [0,0];
L = 100; 
area = L*L;
r0 = 10;

lambda = 0.01; % rwp density


noOfIter = 1000;
distTravelled = zeros(1, noOfIter);
noOfLegs = zeros(1,noOfIter);

for i = 1:noOfIter
    
    currentPosition = origin;
    
    while (norm(currentPosition-origin)<r0)&&(currentPosition(2)>=0)

        N = poissrnd(lambda*area); % no. of AUs
        p = unifrnd(-L/2,L/2,N,2);
        noOfLegs(i) = noOfLegs(i)+1;

        distances = sqrt(sum((p-repmat(currentPosition,length(p),1))'.^2));
        minDist = min(distances);
        nextPosIndex = find(distances == minDist);
        %currentPosition
        nextPosition = p(nextPosIndex,:);

        currentPosition = nextPosition;
        distTravelled(i) = distTravelled(i) + minDist;

    end
    
end

averageDist = sum(distTravelled)/noOfIter
