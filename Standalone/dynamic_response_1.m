function f=dynamic_response_1(fun,N,M,V,first_F1,second_F1,last_pop,t,min_range,max_range,F1_last_Gmax)   % 均为字典
empty_individual.position=[];     
empty_individual.cost=[];        
empty_individual.Rank=[];
empty_individual.DominationSet=[];
empty_individual.DominatedCount=[];
empty_individual.Fitness=[];

 last_pop_fea=[];
last_pop_cv=CV(fun,last_pop,t);
for i=1:length(last_pop_cv)
    if last_pop_cv(i)==0
        last_pop_fea=[last_pop_fea;last_pop(i,:)];
    end
end
if size(last_pop_fea,1)>N/4    
    last_pop_pf=initialize_pf(fun,N, M, V,last_pop_fea,t);
    pop1=repmat(empty_individual,size(last_pop_fea,1),1);   
    for i=1:size(last_pop_fea,1)
        pop1(i).position=last_pop_fea(i,:);
         pop1(i).cost=last_pop_pf(:,i);
    end
    [pop1, F]=NonDominatedSorting(pop1);
    pop1=CalcCrowdingDistance(pop1,F);
    pop1=SortPopulation_2(pop1);

    P1=[];
    for i=1:N/4
     P1=[P1;pop1(i).position];
    end
else      
    P1=[];
    P1=[P1;last_pop_fea];
    flag=find(last_pop_cv==0);
    last_pop(flag,:)=[];            
    last_pop_cv(last_pop_cv==0)=[];
    [~,index]=sort(last_pop_cv,'ascend');
    for i=1:N/4-size(last_pop_fea,1)
        P1=[P1;last_pop(index(i),:)];
    end
end
    

    P2=initialize_variables(fun,N/4, M, V, min_range,max_range,t);  %初始种群
  
    first_position=[];
    second_position=[];
    for i=1:size(first_F1,1)
      first_position=[first_position;first_F1(i).position];     % t-1
      second_position=[second_position;second_F1(i).position];   % t
    end
    center1=sum(first_position)/size(first_position,1);
    center2=sum(second_position)/size(second_position,1);
    direction_vector=center2-center1;
    %%%%计算方差%%%%
    deviation=0;
    for i=1:M
        for k=1:size(F1_last_Gmax,1)
            F_g1=[second_F1.cost];
            F_g=second_F1(k).cost;
            F_last_g=F1_last_Gmax(k).cost;
          individual_deviation=(F_g(i)-F_last_g(i))/(max(F_g1(i,:))-min(F_g1(i,:)));
          miu=abs(individual_deviation)/size(F1_last_Gmax,1);
           deviation= deviation+abs(individual_deviation-miu);
        end
    end
    variance=deviation/(M*size(F1_last_Gmax,1));

    P3=[];
 

for i=1:size(second_position,1)
    x=second_position(i,:);
    noise=variance*randn(1,V);
    new_x=x+direction_vector+noise;
     for jj=1:V
               if(new_x(jj)>max_range(jj))
                    new_x(jj)=max_range(jj);
               elseif(new_x(jj)<min_range(jj))
              new_x(jj)=min_range(jj);
               end
      end
    P3=[P3;new_x];
end
     f=[P1;P2;P3];   % 个数是N个 
end
    
    
