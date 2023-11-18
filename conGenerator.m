function [offX] = conGenerator(popX,popC,minVar,maxVar)
    [popsize,n] = size(popX);
    C = size(popC,2);
    popCsum = sum(max(popC,0),2);
    feaInd = find(popCsum==0);
    offX = popX;
    for i = 1:popsize
        if popCsum(i,1)>0
            l = rand;
            if l <= 1/3
                F  = .6;
            elseif l <= 2/3
                F= 0.8;
            else
                F = 1.0;
            end
            l=rand;
            if l <= 1/3
                CR  = .1;
            elseif l <= 2/3
                CR = 0.2;
            else
                CR = 1.0;
            end
            if ~isempty(feaInd)
                ranFeaInd = unidrnd(length(feaInd));
                bestInd = feaInd(ranFeaInd);
            else
                betterInd = find(popCsum<popCsum(i,1));
                if ~isempty(betterInd)
                    ranBetterInd = unidrnd(length(betterInd));
                    bestInd = betterInd(ranBetterInd);
                else
                    bestInd = unidrnd(popsize-1)+1;
                end
    %             [~,sortInd] = sort(popCsum);
    %             if sortInd(1)==i
    %                 bestInd = sortInd(2);
    %             else
    %                 bestInd = sortInd(1);
    %             end
            end
            % DE/current-to-ranFea/1
            r1 = unidrnd(popsize);
            while (r1==i) || (r1 == bestInd)
                r1 = unidrnd(popsize);
            end
            r2 = unidrnd(popsize);
            while (r2==i) || (r2==r1) || (r2==bestInd)
                r2 = unidrnd(popsize);
            end
            v = popX(i,:)+F*(popX(bestInd,:)-popX(i,:))+F*(popX(r1,:)-popX(r2,:));
            w = find(v<minVar(1,:));
            if ~isempty(w)
                l = rand;
                if l < (1-length(feaInd)/popsize)^3
                    v(1,w) = 2*minVar(1,w)-v(1,w);
                    w1 = find(v(1,w)>maxVar(1,w));
                    if ~isempty(w1)
                        v(1,w(w1)) = minVar(1,w);
                    end
                else
                    v(1,w) = maxVar(1,w);
                end
            end
            y = find(v>maxVar(1,:));
            if ~isempty(y)
                l = rand;
                if l < (1-length(feaInd)/popsize)^3
                    v(1,y) = 2*maxVar(1,y)-v(1,y);
                    y1 = find(v(1,y)<minVar(1,y));
                    if ~isempty(y1)
                        v(1,y(y1)) = minVar(1,y);
                    end
                else
                    v(1,y) = maxVar(1,y);
                end
            end
            t = rand(1,n)<CR;
            j_rand = floor(rand*n)+1;
            t(1,j_rand)=1;
            t_ = 1-t;
            offX(i,:) = t.*v+t_.*popX(i,:);
        end
    end
end

