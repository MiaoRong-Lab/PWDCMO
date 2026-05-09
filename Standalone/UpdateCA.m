function f=UpdateCA(CA,Q,V,M,W)

empty_individual.position=[];     
empty_individual.cost=[];        
empty_individual.Rank=[];
empty_individual.DominationSet=[];
empty_individual.DominatedCount=[];
empty_individual.Fitness=[];


S=[];       
Sc=[];     
SI=[];      
Hc=[CA;Q];
N=size(W,1);
for i=1:size(Hc,1)
    if Hc(i,end)==0
        Sc=[Sc;Hc(i,:)];
    else
        SI=[SI;Hc(i,:)];
    end
end

if size(Sc,1)==N 
    f=Sc(:,1:V);    
elseif size(Sc,1)>N   
    fea=repmat(empty_individual,size(Sc,1),1);   
    for i=1:size(Sc,1)
        fea(i).position=Sc(i,1:V);
        fea(i).cost=Sc(i,V+1:V+M)';
    end
    % 非支配排序
[fea, F]=NonDominatedSorting(fea);
% 计算拥挤度
fea=CalcCrowdingDistance(fea,F);
% 对个体排序
[fea, F]=SortPopulation_2(fea);

for i=1:N
S=[S;fea(i).position];
end
f=S;
elseif size(Sc,1)<N
     S=[S,Sc];
    f1=SI(:,end);    
    PopObj = SI(:,V+1:V+M);
    PopNorm = sqrt(sum(PopObj.^2,2));
    WNorm = sqrt(sum(W.^2,2))';
    [~,Region_SI] = max(PopObj*W'./max(PopNorm*WNorm,eps),[],2);
    Z = min(SI(:,V+1:V+M),[],1) ;
    f2=max(abs(SI(:,V+1:V+M)-repmat(Z,size(SI,1),1))./W(Region_SI,:),[],2);
    PopObj=[f1 f2];    
    infea=repmat(empty_individual,size(SI,1),1);   
    for i=1:size(SI,1)
        infea(i).position=SI(i,1:V);
        infea(i).cost=PopObj(i,:)';
    end
    % 非支配排序
[infea, F]=NonDominatedSorting(infea);
% 计算拥挤度
infea=CalcCrowdingDistance(infea,F);
% 对个体排序
[infea, F]=SortPopulation_2(infea);
rank=1;
F_=infea(F{rank}); 
num=size(Sc,1)+size(F_,1);
while num<N
    rank=rank+1;    
F_=infea(F{rank});   
num=num+size(F_,1);
end
inS1=infea(1:num-size(Sc,1)-size(F_,1)); 
S=S(:,1:V); 
for i=1:size(inS1,1)
    S=[S;inS1(i).position];      
end

    delete_n=num-N;    
    if delete_n==0
        for i=1:size(F_,1)
    S=[S;F_(i).position];  
        end
        f=S;
    else
   objective=[F_.cost];     
    cv=objective(1,:);  
[~,index]=sort(cv,'descend');
index=index(1:delete_n);
F_last=[];    
for i=1:size(F_,1)
    F_last=[F_last;F_(i).position];  
end
for i=1:numel(index)
    a=index(i);     
    F_last(a,:)=0;    
    F_last(sum(F_last,2)==0,:)=[];  
end
 S=[S;F_last];   
    f=S;    
    end
end
end
