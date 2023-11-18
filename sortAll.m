function [rankX,rankY,rankC] = sortAll(popX,popY,popC)
    NP = size(popX,1);
    D = size(popX,2);
    C = size(popC,2);
    [feaInd,infeaInd] = judgeFeasible(popC);
    if length(feaInd) ~= 0
        [sortedFea,sortedFeaInd] = sort(popY(feaInd));
    end
    if length(infeaInd) ~= 0
        popCinfeaSum = sum(max(popC(infeaInd,:),0),2);
        [sortedInfea,sortedInfeaInd] = sort(popCinfeaSum);
    end
    if length(feaInd) ~= 0 && length(infeaInd) ~= 0
        rankX = [popX(feaInd(sortedFeaInd),:);popX(infeaInd(sortedInfeaInd),:)];
        rankY = [popY(feaInd(sortedFeaInd),:);popY(infeaInd(sortedInfeaInd),:)];
        rankC = [popC(feaInd(sortedFeaInd),:);popC(infeaInd(sortedInfeaInd),:)];
%         rankV = [popV(feaInd(sortedFeaInd),:);popV(infeaInd(sortedInfeaInd),:)];
    elseif length(feaInd) ~= 0 && length(infeaInd) == 0
        rankX = popX(feaInd(sortedFeaInd),:);
        rankY = popY(feaInd(sortedFeaInd),:);
        rankC = popC(feaInd(sortedFeaInd),:);
%         rankV = popV(feaInd(sortedFeaInd),:);
    elseif length(feaInd) == 0 && length(infeaInd) ~= 0
        rankX = popX(infeaInd(sortedInfeaInd),:);
        rankY = popY(infeaInd(sortedInfeaInd),:);
        rankC = popC(infeaInd(sortedInfeaInd),:);
%         rankV = popV(infeaInd(sortedInfeaInd),:);
    end
end