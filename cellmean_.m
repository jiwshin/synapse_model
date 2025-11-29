function [outw]  = cellmean_(inw)
    Ncell = length(inw);
    summ = inw{1};
    for i = 2:Ncell
        summ = summ + inw{i};
    end
    outw = summ/Ncell;
end