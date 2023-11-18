function [objF,conV] = fitness(problemSetNum, x, problem, aaa)
    if problemSetNum == 2010
        [objF,conV]=fitness2010(x,problem);
    elseif problemSetNum == 2006
        [objF,conV]=fitness2006(x,problem,aaa);
    end
end

