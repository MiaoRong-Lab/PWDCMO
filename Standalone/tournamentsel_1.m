function p=tournamentsel_1(pop)

n=numel(pop);
s=randperm(n,2);     % 随机选择两个个体
p1=pop(s(1));
p2=pop(s(2));
if p1.Rank<p2.Rank
    p=p1;
elseif p1.Rank>p2.Rank
        p=p2;
elseif p1.Rank==p2.Rank
    if p1.Fitness>p2.Fitness
        p=p1;
    elseif p1.Fitness<p2.Fitness
        p=p2;
    elseif p1.Fitness== p2.Fitness
        s1=s(randperm(length(s),1));
        p=pop(s1);

    end
end
end
