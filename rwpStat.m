origin = [0,0];
L = 1000; 
area = L*L;
lambdaS = 0.001:0.0099:0.1; % rwp density

meanLegLengthSim = zeros(1,numel(lambdaS));
meanLegLength = zeros(1,numel(lambdaS));
noOfLegs = 100;

noOfPaths = 10;
totalDist = zeros(1,noOfPaths);

for k = 1:numel(lambdaS)
    lambda = lambdaS(k);

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
meanLegLengthSim(k) = sum(totalDist)/(noOfLegs*noOfPaths);
meanLegLength(k) = 1/(2*sqrt(lambda));
end

plot(lambdaS,meanLegLength,'r',lambdaS,meanLegLengthSim);
legend('analytic','simulation');
title('Mean Transition Length');
xlabel('\lambda');
ylabel('E[L]');