function [NewSequence1 , NewSequence2, NewSequence3, NewSequence4] = Forage_Best (OldSequence, nghk)
    nghk = ceil(nghk);
            NewSequence1 = Swap (OldSequence, nghk);
            NewSequence2 = Reversion (OldSequence, nghk);
            NewSequence3 = Insertion (OldSequence, nghk);
            NewSequence4 = Multi_Insert (OldSequence, nghk);
end

function NewSequence = Swap (OldSequence, nghk)
    k = numel(OldSequence);
    i1 = randi([1 k]);
    i2 = i1 + randi([1 nghk]);
    i2(i2>k) = k;
    NewSequence = OldSequence;
    NewSequence([i1 i2]) = OldSequence([i2 i1]);
end

function NewSequence = Reversion (OldSequence, nghk)
    k = numel (OldSequence);
    i1 = randi ([1 k]);
    i2 = i1 + randi ([1 nghk]);
    i2(i2>k) = k;

    NewSequence = OldSequence;
    NewSequence(i1:i2) = OldSequence(i2:-1:i1);
end

function NewSequence = Insertion (OldSequence, nghk)
    k = numel (OldSequence);
    a = randi(2);
    switch a
        case 1
            i1 = randi ([1 k]);
        case 2
            i1 = randi ([1 k-1]);
            i1 = [i1 i1+1];
    end
    i2 = i1(end) + randi([-nghk nghk]);
    i2(i2>k) = k;
    i2(i2<1) = 1;

    if i1<i2
        NewSequence = [OldSequence(1:i1-1) OldSequence(i1+1:i2) OldSequence([i1]) OldSequence(i2+1:end)];
    else
        NewSequence = [OldSequence(1:i2) OldSequence([i1]) OldSequence(i2+1:i1-1) OldSequence(i1+1:end)];
    end
end

function NewSequence = Multi_Insert (OldSequence, nghk)
    nmc = randi([2 5],1);
    nc = numel(OldSequence);
    m1 = randi([1 (nc-(nmc-1))],1); %Generate rnadom m1
    
    for i=1:nmc
        mc1(i)=(m1+i-1);
    end
    %Select a random area m2 inthe ngh excep for mc1
    m2_1 = randi([1 m1]);
    m2_2 = randi([m1 nc]);
switch randi([1 2])
    case 1
        m2 = m2_1;
    case 2
        m2 = m2_2;
end
    %Generage a random number t (1-2)
    t = randi ([1 2],1);
    if t==1
        if m1<m2
        NewSequence = OldSequence([1:m1-1   m1+nmc:m2   m1:m1+nmc-1 m2+1:nc]);
        else
        NewSequence = OldSequence([1:m2   m1:m1+nmc-1   m2+1:m1-1 m1+nmc:nc]);
        end
    else
        if m1<m2
            NewSequence = OldSequence([1:m1-1   m1+nmc:m2   m1+nmc-1:-1:m1 m2+1:nc]);
        else
            NewSequence = OldSequence([1:m2   m1+nmc-1:-1:m1   m2+1:m1-1 m1+nmc:nc]);
        end
    end
end