function [popY_pen] = calPenaltyY(popY,popC)
    [popsize,C] = size(popC);
    [feaInd,infeaInd] = judgeFeasible(popC);
    popsize=length(feaInd)+length(infeaInd);
    feaRatio = length(feaInd)/popsize;
    popY_norm = (popY-min(popY))./(max(popY)-min(popY));
    popC_max = max(popC,[],1);
    popC_norm = ones(popsize,1)*1e16;
    for i = 1:popsize
        popC_norm(i,1) = 0;
        for j = 1:C
            if popC_max(j)>0
                popC_norm(i,1) = popC_norm(i,1)+popC(i,j)/popC_max(j);
            end
        end
        popC_norm(i,1) = popC_norm(i,1)/C;
    end
    if feaRatio == 0
        dx = popC_norm;
        Xx = zeros(popsize,1);
    else
        dx = sqrt(popC_norm.*popC_norm+popY_norm.*popY_norm);
        Xx = popC_norm;
    end
    Yx = ones(popsize,1)*1e16;
    for ii = 1:popsize
        if ismember(ii,feaInd)
            Yx(ii,:) = 0;
        else
            Yx(ii,:) = popY_norm(ii,:);
        end
    end
    px = (1-feaRatio).*Xx+feaRatio.*Yx;
    popY_pen = dx+px;
end

