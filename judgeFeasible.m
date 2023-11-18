function [feaInd,infeaInd] = judgeFeasible(originC)
    N = size(originC,1);
    C = size(originC,2);
    feaInd = [];
    infeaInd = [];
    for i = 1:N
        isFeasible = 1;
        for j = 1:C
            if originC(i,j) > 0
                isFeasible = 0;
                break;
            end
        end
        if isFeasible == 1
            feaInd = [feaInd;i];
        else
            infeaInd = [infeaInd;i];
        end
    end
end

