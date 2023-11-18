function [rectifiedY,rectifiedC] = neighbor(offConX,offConY,offConC,archiveX,archiveY,archiveCsum,modelCon,modelObj,flagCon)
    omega = 2; 
    [popsize,n] = size(offConX);
    %%%%%%%neighborhood information
    for i = 1:popsize
        neighDis = zeros(size(archiveX,1),1);
        minNeighDis = 1e16;
        minNeighDisInd = 1e4;
        for j = 1:size(archiveX,1)
            for k = 1:n
                neighDis(j,1) = neighDis(j,1)+(offConX(i,k)-archiveX(j,k))*(offConX(i,k)-archiveX(j,k));
            end
            neighDis(j,1) = sqrt(neighDis(j,1));
            if neighDis(j,1)<minNeighDis
                minNeighDis = neighDis(j,1);
                minNeighDisInd = j;
            end
        end
%         [minNeighDis,minNeighDisInd] = min(neighDis);
        if minNeighDis~=1e16 && minNeighDisInd~=1e4
            [mmeanY,~,mmseY] = predictor(archiveX(minNeighDisInd,:),modelObj);
            testErrorY = archiveY(minNeighDisInd,:)-(mmeanY-omega.*mmseY);
            offConY(i,1) = offConY(i,1)+testErrorY;
            if flagCon==1
                [mmeanC,~,mmseC] = predictor(archiveX(minNeighDisInd,:),modelCon);
                testErrorC = archiveCsum(minNeighDisInd,:)-(mmeanC-omega.*mmseC);
                offConC(i,1) = offConC(i,1)+testErrorC;
            end
        end
    end
    %%%%%%%neighborhood information
    rectifiedY = offConY;
    rectifiedC = offConC;
end

