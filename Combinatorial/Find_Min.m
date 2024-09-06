function [BestP, BestC] = Find_Min(Pos1, Pos2, Pos3, Pos4, Cost1, Cost2, Cost3, Cost4)
    if Cost1<Cost2
        Best1P = Pos1;
        Best1C = Cost1;
    else
        Best1P = Pos2;
        Best1C = Cost2;
    end
    if Cost3<Cost4
        Best2P = Pos3;
        Best2C = Cost3;
    else
        Best2P = Pos4;
        Best2C = Cost4;
    end
    if Best1C<Best2C
        BestP = Best1P;
        BestC = Best1C;
    else
        BestP = Best2P;
        BestC = Best2C;
    end
end