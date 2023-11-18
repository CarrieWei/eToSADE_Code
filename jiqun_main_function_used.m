function [evolveSolution,evolveConstrain,run_time,feasiRatio] = jiqun_main_function_used(problemSetNum,problem,time,popsize,n,totalFES,totalTime,minVar,maxVar,aaa)
    format long g;
    format compact;
    addpath(genpath('./dace/'))
    set(0,'RecursionLimit',2000000);
    'SADE'
    filename = ['process/' sprintf('%d_%d_n%d_NP%d_FEs%d_runs%d_%d_process.txt', problemSetNum, problem,n,popsize,totalFES, totalTime,time)];
    fid = fopen(filename,'w');
    filename_gen = ['process/' sprintf('conflag_%d_%d_n%d_NP%d_FEs%d_runs%d_%d.csv', problemSetNum, problem,n,popsize,totalFES, totalTime,time)];
    fid_gen = fopen(filename_gen,'w');
    numTrain = 500;
    classNum = round(numTrain / 3);
    numEvalu = 2;
    omega = 2;
    Nrand = 5;
    threshold = 0;
    maxGen = totalFES/numEvalu;
%     problemSetNum = 2006;
    run_time = zeros(totalTime,1);
    feasiRatio = zeros(totalFES, totalTime);
    evolveSolution = ones(totalFES,totalTime)*1e16;
    evolveConstrain = ones(totalFES,totalTime)*1e16;
    start_time = cputime;
    archiveX = lhsamp(numTrain,n);
    archiveX = archiveX.*repmat(maxVar-minVar,numTrain,1)+repmat(minVar,numTrain,1);
    [archiveY,archiveC] = fitness(problemSetNum,archiveX,problem,aaa);
    [archiveX,archiveY,archiveC] = sortAll(archiveX,archiveY,archiveC);
    archiveCsum  = sum(max(archiveC,0), 2);
    [archiveFeaInd,archiveInfeaInd] = judgeFeasible(archiveCsum);
    theta=ones(1,n).*0.5;
    modelObj = dacefit(archiveX,archiveY,@regpoly2,@correxp,theta,zeros(1,n)+1e-10,ones(1,n));
    fprintf('Fea %d, Infea %d\n', length(archiveFeaInd),length(archiveInfeaInd));
    fprintf(fid,'%d, %d\n', length(archiveFeaInd),length(archiveInfeaInd));
    flagCon = -1;
    if length(archiveFeaInd) < classNum%isempty(archiveFeaInd) %
        flagCon = 1;
    elseif ~isempty(archiveInfeaInd) && length(archiveFeaInd) >= classNum %~isempty(archiveFeaInd)%
        flagCon = 2;
    elseif isempty(archiveInfeaInd)
        flagCon = 3;
    end
    if flagCon == 1 
        modelCon = dacefit(archiveX,archiveCsum,@regpoly2,@correxp,theta,zeros(1,n)+1e-10,ones(1,n));
    elseif flagCon == 2
        trainC = archiveCsum;
        trainC(archiveFeaInd,:) = -1;
        trainC(archiveInfeaInd,:) = 1;
        modelCon = fitcsvm(archiveX,trainC,'Standardize',true,'KernelFunction','RBF','KernelScale','auto');
        CVSVMModel = crossval(modelCon);
        classLoss = kfoldLoss(CVSVMModel);
        fprintf('classLoss: %f\n', classLoss);
        fprintf(fid,'classLoss: %f\n', classLoss);
    elseif flagCon == 3
        modelCon = [];
    end
    popX = archiveX(1:popsize,:);
    popY = archiveY(1:popsize,:);
    popC = archiveC(1:popsize,:);
    popCsum = sum(max(popC,0), 2);
    FES = 1;
    gen = 0;
    fprintf('run, %d %d; best, %f %f\n', time, FES, popCsum(1,:), popY(1,:));
    fprintf(fid,'%d %d; %f %f\n', time, FES, popCsum(1,:), popY(1,:));
    while FES < totalFES %FES<1%
        fprintf(fid_gen,'%d,%d\n',FES,flagCon);
        %%% Sample Generation
        if problem==6||problem==11
            if FES==1-1||FES==250-1||FES==500-1||FES==750-1||FES==1000-1
                if problem==6
                    if FES==0
                        lbound=[10,0];
                        ubound=[100,100];
    %                 elseif FES==250
    %                     lbound=[10,0];
    %                     ubound=[50,50];
                    else
                        lbound=[13,0];
                        ubound=[16,12];
                    end
                else
                    lbound = minVar;
                    ubound = maxVar;
                end
                FileName = [sprintf('samples/') sprintf('g%d',problem) '/' sprintf('pop_%d_%d_n%d_FES%d_run%d_%d.csv', problemSetNum, problem,n,FES+1,totalTime,time)];
                dlmwrite(FileName, [popX,popY,popC], 'precision', '%.6f');
                NI=100;
                step1 = (ubound(1,1)-lbound(1,1))/(NI-1);
                step2 = (ubound(1,2)-lbound(1,2))/(NI-1);
                sampleX1 = zeros(NI,1);
                sampleX2 = zeros(NI,1);
                for i = 1:NI
                    sampleX1(i,1) = lbound(1,1)+step1*(i-1);
                    sampleX2(i,1) = lbound(1,2)+step2*(i-1);
                end
                resultsC = zeros(NI,NI);
                if flagCon==1
                    for i = 1:NI
                        for j = 1:NI
                            sampleX = [sampleX1(i,1),sampleX2(j,1)];
                            [meanC,~,mseC] = predictor(sampleX,modelCon);
                            resultsC(j,i) = meanC-omega.*mseC;
                        end
                    end
                    FileName = [sprintf('samples/') sprintf('g%d',problem) '/' sprintf('Stage1_%d_%d_n%d_FES%d_run%d_%d.csv', problemSetNum, problem,n,FES+1,totalTime,time)];
                    dlmwrite(FileName, resultsC, 'precision', '%.6f');
                elseif flagCon==2
                    for i = 1:NI
                        for j = 1:NI
                            sampleX = [sampleX1(i,1),sampleX2(j,1)];
                            resultsC(j,i) = predict(modelCon,sampleX);
                        end
                    end
                    FileName = [sprintf('samples/') sprintf('g%d',problem) '/' sprintf('Stage2_%d_%d_n%d_FES%d_run%d_%d.csv', problemSetNum, problem,n,FES+1,totalTime,time)];
                    dlmwrite(FileName, resultsC, 'precision', '%.6f');
                end
            end
        end
        %%% Sample Generation
        % the first operator
        feaInd = find(archiveCsum<=0);
        if length(feaInd) < classNum
            fprintf('1.Feasible Exploration\n');
            offConX = conGenerator(popX,popC,minVar,maxVar);
            meanY = [];
            mseY = [];
            meanC = [];
            mseC = [];
            if size(offConX,1)==1
                [meanY,~,mseY] = predictor(offConX,modelObj);
            else
                [meanY,mseY] = predictor(offConX,modelObj);
            end
            offConY = meanY-omega.*mseY;
            if flagCon == 1 
                if size(offConX,1)==1
                    [meanC,~,mseC] = predictor(offConX,modelCon);
                else
                    [meanC,mseC] = predictor(offConX,modelCon);
                end
                offConC = meanC-omega.*mseC;
            elseif flagCon == 2
                offConC = predict(modelCon,offConX);
            elseif flagCon == 3
                offConC = zeros(size(offConX,1),1);
            end 
            [offConY,offConC] = neighbor(offConX,offConY,offConC,archiveX,archiveY,archiveCsum,modelCon,modelObj,flagCon);
            offConC = max(offConC,0);
            offConFeaInd = find(offConC<=0); 
            count = 0;
            while length(offConFeaInd)<popsize
                if count == 100
                    break;
                end
                count = count+1;
                offConX2 = conGenerator(offConX,offConC,minVar,maxVar);
                meanY2 = [];
                mseY2 = [];
                meanC2 = [];
                mseC2 = [];
                if size(offConX2,1)==1
                    [meanY2,~,mseY2] = predictor(offConX2,modelObj);
                else
                    [meanY2,mseY2] = predictor(offConX2,modelObj);
                end
                offConY2 = meanY2-omega.*mseY2;
                if flagCon == 1 
                    if size(offConX2,1)==1
                        [meanC2,~,mseC2] = predictor(offConX2,modelCon);
                    else
                        [meanC2,mseC2] = predictor(offConX2,modelCon);
                    end
                    offConC2 = meanC2-omega.*mseC2;
                elseif flagCon == 2
                    offConC2 = predict(modelCon,offConX2);
                elseif flagCon == 3
                    offConC2 = zeros(size(offConX2,1),1);
                end
                [offConY2,offConC2] = neighbor(offConX2,offConY2,offConC2,archiveX,archiveY,archiveCsum,modelCon,modelObj,flagCon);
                offConC2 = max(offConC2,0);
                %%% Crowding-based infeasible replacement
                for i = 1:popsize
                    minDis = 1e16;
                    minDisInd = 0;
                    offConInfeaInd = find(offConC~=0);
                    if ~isempty(offConInfeaInd)
                        for iii = 1:length(offConInfeaInd)
                            ii = offConInfeaInd(iii);
                            distance = 0;
                            for j = 1:n
                                distance = distance+(offConX2(i,j)-offConX(ii,j))*(offConX2(i,j)-offConX(ii,j));
                            end
                            distance = sqrt(distance);
                            if distance < minDis
                                minDis = distance;
                                minDisInd = ii;
                            end
                        end
                    else
                        for ii = 1:popsize
                            distance = 0;
                            for j = 1:n
                                distance = distance+(offConX2(i,j)-offConX(ii,j))*(offConX2(i,j)-offConX(ii,j));
                            end
                            distance = sqrt(distance);
                            if distance < minDis
                                minDis = distance;
                                minDisInd = ii;
                            end
                        end
                    end
                    offConX(minDisInd,:) = offConX2(i,:);
                    offConY(minDisInd,:) = offConY2(i,:);
                    offConC(minDisInd,:) = offConC2(i,:);
                end
                offConFeaInd = find(offConC<=0);    
            end
            fprintf('count1 %d\n', count);
            fprintf(fid,'%d\n', count);
            if isempty(offConFeaInd)
                [~,selectInd] = min(offConC);
            else
                [~,worstInd] = min(offConY(offConFeaInd,:));
                selectInd = offConFeaInd(worstInd);
            end
            selectX = offConX(selectInd,:);
            if ~ismember(selectX,archiveX,'rows')
                [selectY,selectC] = fitness(problemSetNum,selectX,problem,aaa);
                FES = FES+1;
                fprintf('Selecting offConX\n');
                archiveX = [archiveX;selectX];
                archiveY = [archiveY;selectY];
                archiveC = [archiveC;selectC];
                archiveCsum = sum(max(archiveC,0), 2);
            end
            % the evolution of objective
            offConY_pen = calPenaltyY(offConY,offConC);
            offObjX = objGenerator(offConX,offConY_pen,minVar,maxVar);
        else
            popY_pen = calPenaltyY(popY,popC);
            offObjX = objGenerator(popX,popY_pen,minVar,maxVar);
        end
        meanY = [];
        mseY = [];
        meanC = [];
        mseC = [];
        if size(offObjX,1)==1
            [meanY,~,mseY] = predictor(offObjX,modelObj);
        else
            [meanY,mseY] = predictor(offObjX,modelObj);
        end
        offObjY = meanY-omega.*mseY;
        if flagCon == 1 
            if size(offObjX,1)==1
                [meanC,~,mseC] = predictor(offObjX,modelCon);
            else
                [meanC,mseC] = predictor(offObjX,modelCon);
            end
            offObjC = meanC-omega.*mseC;
        elseif flagCon == 2
            offObjC = predict(modelCon,offObjX);
        elseif flagCon == 3
            offObjC = zeros(size(offObjX,1),1);
        end
        [offObjY,offObjC] = neighbor(offObjX,offObjY,offObjC,archiveX,archiveY,archiveCsum,modelCon,modelObj,flagCon);
        offObjC = max(offObjC,0);
        offObjFeaInd = find(offObjC<=0);
        offObjY_pen = calPenaltyY(offObjY,offObjC);
        count = 0;
        while length(offObjFeaInd) <= 2
            if count == 100
                break;
            end
            count = count+1;
            offObjX2 = objGenerator(offObjX,offObjY_pen,minVar,maxVar);
            meanY2 = [];
            mseY2 = [];
            meanC2 = [];
            mseC2 = [];
            if size(offObjX2,1)==1
                [meanY2,~,mseY2] = predictor(offObjX2,modelObj);
            else
                [meanY2,mseY2] = predictor(offObjX2,modelObj);
            end
            offObjY2 = meanY2-omega.*mseY2;
            if flagCon == 1 
                if size(offObjX2,1)==1
                    [meanC2,~,mseC2] = predictor(offObjX2,modelCon);
                else
                    [meanC2,mseC2] = predictor(offObjX2,modelCon);
                end
                offObjC2 = meanC2-omega.*mseC2;
            elseif flagCon == 2
                offObjC2 = predict(modelCon,offObjX2);
            elseif flagCon == 3
                offObjC2 = zeros(size(offObjX2,1),1);
            end
            [offObjY2,offObjC2] = neighbor(offObjX2,offObjY2,offObjC2,archiveX,archiveY,archiveCsum,modelCon,modelObj,flagCon);
            offObjC2 = max(offObjC2,0);
            offObjY2_pen = calPenaltyY(offObjY2,offObjC2);
            for i = 1:popsize
                if offObjY2_pen(i,:)<offObjY_pen(i,:)
                    offObjX(i,:) = offObjX2(i,:);
                    offObjY(i,:) = offObjY2(i,:);
                    offObjC(i,:) = offObjC2(i,:);
                end
            end
			offObjY_pen = calPenaltyY(offObjY,offObjC);
            offObjFeaInd = find(offObjC<=0);
        end
        fprintf('count2 %d\n', count);
        fprintf(fid,'%d\n', count);
        if flagCon==1||flagCon==2
            offObjY_pen = calPenaltyY(offObjY,offObjC);
            [~,sortOffObjInd] = sort(offObjY_pen);
        else
            [~,sortOffObjInd] = sort(offObjY);
        end
        selectInd = sortOffObjInd(1);
        selectX = offObjX(selectInd,:);
        if length(feaInd) >= classNum
            offObjX(selectInd,:) = [];
            offObjY(selectInd,:) = [];
            offObjC(selectInd,:) = [];
            offFeaInd = find(offObjC<=0);
            if ~isempty(offFeaInd)
                [~,bestInd] = min(offObjY(offFeaInd));
                selectInd2 = offFeaInd(bestInd);
            else
                [~,selectInd2] = min(offObjC);
            end
            selectX = [selectX;offObjX(selectInd2,:)];
        end
        for dd = 1:length(selectInd)
            if ~ismember(selectX(dd,:),archiveX,'rows')
                [selectY(dd,:),selectC(dd,:)] = fitness(problemSetNum,selectX(dd,:),problem,aaa);
                FES = FES+1;
                fprintf('Selecting offObjX\n');
                archiveX = [archiveX;selectX(dd,:)];
                archiveY = [archiveY;selectY(dd,:)];
                archiveC = [archiveC;selectC(dd,:)];
                archiveCsum = sum(max(archiveC,0), 2);
            end
        end
        [archiveX,archiveY,archiveC] = sortAll(archiveX,archiveY,archiveC);
        archiveCsum  = sum(max(archiveC,0), 2);
        modelObj = dacefit(archiveX(1:numTrain,:),archiveY(1:numTrain,:),@regpoly2,@correxp,theta,zeros(1,n)+1e-10,ones(1,n));
        [archiveFeaInd,archiveInfeaInd] = judgeFeasible(archiveCsum);
        fprintf('Fea %d, Infea %d\n', length(archiveFeaInd),length(archiveInfeaInd));
        fprintf(fid,'%d, %d\n', length(archiveFeaInd),length(archiveInfeaInd));
        flagCon = -1;
        if length(archiveFeaInd) < classNum
            flagCon = 1;
        elseif ~isempty(archiveInfeaInd) && length(archiveFeaInd) >= classNum
            flagCon = 2;
        elseif isempty(archiveInfeaInd)
            flagCon = 3;
        end
        % boundary training data set selection
        if ~isempty(archiveInfeaInd)
            feaCanX = archiveX(archiveFeaInd,:);
            infeaCanX = archiveX(archiveInfeaInd,:);
            if isempty(archiveFeaInd) %全是不可行解，选做好的300个个体来做回归
                conTrainX = archiveX(1:numTrain,:);
                conTrainCsum = archiveCsum(1:numTrain,:);
            elseif ~isempty(archiveFeaInd) && length(archiveFeaInd)<=numTrain-classNum %存在可行解但可行解个数少于做分类的个数，选择所有的可行解和离可行解最近的不可行解做回归
                Distances = zeros(length(archiveFeaInd),length(archiveInfeaInd));
                for d1 = 1:length(archiveFeaInd)
                    for d2 = 1:length(archiveInfeaInd)
                        for d3 = 1:n
                            Distances(d1,d2) = Distances(d1,d2)+(feaCanX(d1,d3)-infeaCanX(d2,d3))*(feaCanX(d1,d3)-infeaCanX(d2,d3));
                        end
                        Distances(d1,d2) = sqrt(Distances(d1,d2));
                    end
                end
                mergeDisVec = max(Distances,[],1);%为防止扎堆
                [~,sortVecInd] = sort(mergeDisVec);
                numInfea = numTrain - length(archiveFeaInd);
                selInfeaInd = archiveInfeaInd(sortVecInd(1:numInfea));
                conTrainX = [archiveX(archiveFeaInd,:);archiveX(selInfeaInd,:)];
                conTrainCsum = [archiveCsum(archiveFeaInd,:);archiveCsum(selInfeaInd,:)];
            elseif ~isempty(archiveFeaInd) && length(archiveInfeaInd)<=classNum %存在可行解但不可行解个数少于做分类的个数，选择所有的不可行解和离不可行解最近的可行解做回归
                Distances = zeros(length(archiveInfeaInd),length(archiveFeaInd));
                for d1 = 1:length(archiveInfeaInd)
                    for d2 = 1:length(archiveFeaInd)
                        for d3 = 1:n
                            Distances(d1,d2) = Distances(d1,d2)+(infeaCanX(d1,d3)-feaCanX(d2,d3))*(infeaCanX(d1,d3)-feaCanX(d2,d3));
                        end
                        Distances(d1,d2) = sqrt(Distances(d1,d2));
                    end
                end
                mergeDisVec = max(Distances,[],1);%为防止扎堆
                [~,sortVecInd] = sort(mergeDisVec);
                numFea = numTrain - length(archiveInfeaInd);
                selFeaInd = archiveFeaInd(sortVecInd(1:numFea));
                conTrainX = [archiveX(archiveInfeaInd,:);archiveX(selFeaInd,:)];
                conTrainCsum = [archiveCsum(archiveInfeaInd,:);archiveCsum(selFeaInd,:)];
            elseif length(archiveFeaInd)>numTrain-classNum && length(archiveInfeaInd)>classNum %可行解与不可行解均大于classNum个，不可行解仅保留classNum个
                conTrainFeaX = archiveX(length(archiveFeaInd)-(numTrain-classNum)+1:length(archiveFeaInd),:);
                conTrainFeaCsum = archiveCsum(length(archiveFeaInd)-(numTrain-classNum)+1:length(archiveFeaInd),:);
                Distances = zeros(size(conTrainFeaX,1),length(archiveInfeaInd));
                for d1 = 1:size(conTrainFeaX,1)
                    for d2 = 1:length(archiveInfeaInd)
                        for d3 = 1:n
                            Distances(d1,d2) = Distances(d1,d2)+(conTrainFeaX(d1,d3)-infeaCanX(d2,d3))*(conTrainFeaX(d1,d3)-infeaCanX(d2,d3));
                        end
                        Distances(d1,d2) = sqrt(Distances(d1,d2));
                    end
                end
                mergeDisVec = max(Distances,[],1);%为防止扎堆
                [~,sortVecInd] = sort(mergeDisVec);
                numInfea = numTrain - size(conTrainFeaX,1);
                selInfeaInd = archiveInfeaInd(sortVecInd(1:numInfea));
                conTrainX = [conTrainFeaX;archiveX(selInfeaInd,:)];
                conTrainCsum = [conTrainFeaCsum;archiveCsum(selInfeaInd,:)];
            end
        end
        if flagCon == 1 
            modelCon = dacefit(conTrainX,conTrainCsum,@regpoly2,@correxp,theta,zeros(1,n)+1e-10,ones(1,n));
        elseif flagCon == 2
            trainC = conTrainCsum;
            ind1 = find(conTrainCsum<=0);
            ind2 = find(conTrainCsum>0);
            trainC(ind1,:) = -1;
            trainC(ind2,:) = 1;
            modelCon = fitcsvm(conTrainX,trainC,'Standardize',true,'KernelFunction','RBF','KernelScale','auto');
            CVSVMModel = crossval(modelCon);
            classLoss = kfoldLoss(CVSVMModel);
            fprintf('classLoss: %f\n', classLoss);
            fprintf(fid,'classLoss: %f\n', classLoss);
        elseif flagCon == 3
        end
        popX = archiveX(1:popsize,:);
        popY = archiveY(1:popsize,:);
        popC = archiveC(1:popsize,:);
%         if isempty(archiveFeaInd)
%             popX = archiveX(1:popsize,:);
%             popY = archiveY(1:popsize,:);
%             popC = archiveC(1:popsize,:);
%         else
%             popX = archiveX(1:popsize-Nrand,:);
%             popY = archiveY(1:popsize-Nrand,:);
%             popC = archiveC(1:popsize-Nrand,:);
%             if length(archiveInfeaInd)>=Nrand
%                 randInd = (archiveFeaInd+1:archiveFeaInd+Nrand);
%             else
%                 randInd = unidrnd(size(archiveX,1)-(popsize-Nrand),[Nrand,1]);
%                 randInd = randInd + (popsize-Nrand); 
%             end
%             popX = [popX;archiveX(randInd,:)];
%             popY = [popY;archiveY(randInd,:)];
%             popC = [popC;archiveC(randInd,:)];
%         end
        popCsum = sum(max(popC,0), 2);
        fprintf('run, %d %d; best, %f %f\n', time, FES, popCsum(1,:), popY(1,:));
        fprintf(fid,'%d %d; %f %f\n', time, FES, popCsum(1,:), popY(1,:));
        [feaIndP,infeaIndP] = judgeFeasible(popC);
        if ~isempty(feaIndP)
            feaRatio = length(feaIndP)/(length(feaIndP)+length(infeaIndP));
        else
            feaRatio = 0;
        end
        [feaIndA,~] = judgeFeasible(archiveC);
        if isempty(feaIndA)
            bestSolution = NaN;
            bestConstrain = min(archiveCsum);
        else
            bestSolution = min(archiveY(feaIndA,:));
            bestConstrain = min(archiveCsum);
        end
        gen = gen+1;
        evolveSolution(FES,time) = bestSolution;
        evolveConstrain(FES,time) = bestConstrain;
        feasiRatio(FES,time) = feaRatio;
    end
    run_time(time,1) = cputime-start_time;
    fclose(fid);
    fclose(fid_gen);
end

