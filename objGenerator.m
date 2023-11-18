function [offX] = objGenerator(popX,popY,minVar,maxVar)
    [popsize,n] = size(popX);
    trialX = ones(popsize,n)*1e16;
%     for i = 1:popsize
%         l = rand;
%         if l <= 1/3
%             F  = .6;
%         elseif l <= 2/3
%             F= 0.8;
%         else
%             F = 1.0;
%         end
%         l=rand;
%         if l <= 1/3
%             CR  = .1;
%         elseif l <= 2/3
%             CR = 0.2;
%         else
%             CR = 1.0;
%         end
%         [~,bestInd] = min(popY);
%         % DE/current-to-best/1
%         r1 = unidrnd(popsize);
%         if length(bestInd)>1
%             bestInd = bestInd(1);
%         end
%         while (r1==i) || (r1==bestInd)
%             r1 = unidrnd(popsize);
%         end
%         r2 = unidrnd(popsize);
%         while (r2==i) || (r2==r1) || (r2==bestInd)
%             r2 = unidrnd(popsize);
%         end
%         v = popX(i,:)+F*(popX(bestInd,:)-popX(i,:))+F*(popX(r1,:)-popX(r2,:));
%         for j = 1:n
%             if v(1,j)>maxVar(1,j) || v(1,j)<minVar(1,j)
%                 v(1,j) = rand()*(maxVar(1,j)-minVar(1,j))+minVar(1,j);
%             end
%         end
%         t = rand(1,n)<CR;
%         j_rand = floor(rand*n)+1;
%         t(1,j_rand)=1;
%         t_ = 1-t;
%         trialX(i,:) = t.*v+t_.*popX(i,:);
%     end
    F = 0.8;
    CR = 0.8;
    [~,bestYind] = min(popY);
    for i = 1:popsize
        indexSet = 1:popsize;
        indexSet(i) = [];
       % indexSet(bestYind) = [];
        temp = floor(rand*(popsize-1))+1;
        index1 = indexSet(temp);
        while index1 == bestYind
            temp = floor(rand*(popsize-1))+1;
            index1 = indexSet(temp);
        end
        temp = floor(rand*(popsize-1))+1;
        index2 = indexSet(temp);
        while index2 == bestYind || index2 == index1
            temp = floor(rand*(popsize-1))+1;
            index2 = indexSet(temp);
        end
        v = popX(bestYind,:)+F.*(popX(index1,:)-popX(index2,:));
        for j = 1:n
            if (v(1,j) < minVar(1,j)) || (v(1,j) > maxVar(1,j))
                v(1,j) = rand()*(maxVar(1,j)-minVar(1,j))+minVar(1,j);
            end
        end
        jRand = floor(rand * n) + 1;
        t = rand(1, n) < CR;
        t(1, jRand) = 1;
        t_ = 1 - t;
        u = t .* v + t_ .* popX(i, : );
        offX(i,:) = u;
    end
%     offX = trialX;
end

