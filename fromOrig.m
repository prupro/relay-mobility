origin = [0,0];
L = 1000; 
area = L*L;

lambda = 0.001; % rwp density

radius1 = 50;
radius2 = 1;
center1 = origin;
center2 = [radius1,0];

noOfIter = 1000;
distTravelled = zeros(1, noOfIter);
noOfLegs = zeros(1,noOfIter);

for i = 1:noOfIter
    
    currentPosition = center2;
    
    while (norm(currentPosition-center1)<=radius1)&&(norm(currentPosition-center2)<=radius2)

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
    distTravelled(i) = distTravelled(i) - minDist +  min( incircleLength(prevPosition, currentPosition, radius1, center1),    incircleLength(prevPosition, currentPosition, radius2, center2) );
    
end

averageDist = sum(distTravelled)/noOfIter
