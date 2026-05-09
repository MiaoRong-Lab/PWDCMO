%% 种群排序
function [pop, F]=SortPopulation_2(pop)
% 基于拥挤度排序
[~, CDSO]=sort([pop.Fitness],'descend');  % 'descend'
pop=pop(CDSO);
[~, RSO]=sort([pop.Rank]);
pop=pop(RSO);
Ranks=[pop.Rank];
MaxRank=max(Ranks);
F=cell(MaxRank,1);
for r=1:MaxRank
    F{r}=find(Ranks==r);
end
end