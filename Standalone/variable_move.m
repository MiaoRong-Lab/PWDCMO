function f= variable_move(fun,N,M,V,t,min_range,max_range,pop)
switch fun
    case 1
    % DCMOP1
s=max(3.5-0.14*t,0.7+0.14*t);       %%%连续-离散-连续
m=max(1.43-0.05*t,0.43+0.05*t);
w=2;
    case 2
% DCMOP2
s=max(2.5-0.05*t,1.5+0.05*t);       % 离散-连续-离散
m=max(1.16-0.075*t,-0.43+0.075*t);
w=2;

    case 3
%DCMOP3
s=min(2.1+0.14*t,4.9-0.14*t);       % 离散-连续-离散
m=min(0.93+0.05*t,1.93-0.05*t);
w=2;

 case 4
% %% DCMOP4
s=min(2+0.05*t,3-0.05*t);       % 连续-离散-连续
m=min(0.41+0.075*t,1.91-0.075*t);
w=2;

 case 5
% %% DCMOP5
s=max(3.5-0.14*t,0.7+0.14*t);       %连续-离散-连续
m=max(1.43-0.05*t,0.43+0.05*t);
w=6*sin(0.2*pi*(t+1));

 case 6
% % DCMOP6
s=max(2.5-0.05*t,1.5+0.05*t);       % 离散-连续-离散
m=max(1.16-0.075*t,-0.43+0.075*t);
w=6*sin(0.2*pi*(t+1));

 case 7
% %% DCMOP7
s=min(2.1+0.14*t,4.9-0.14*t);       % 离散-连续-离散
m=min(0.93+0.05*t,1.93-0.05*t);
w=6*sin(0.2*pi*(t+1));

 case 8
% %% DCMOP8
s=min(2+0.05*t,3-0.05*t);       % 连续-离散-连续
m=min(0.41+0.075*t,1.91-0.075*t);
w=6*sin(0.2*pi*(t+1));
end

empty_individual.position=[];     
empty_individual.cost=[];        
empty_individual.Rank=[];
empty_individual.DominationSet=[];
empty_individual.DominatedCount=[];
empty_individual.Fitness=[];


feasible_position=[];
feasible_cost=[];
infeasible_position=[];
final_pop=[];
   for l=1:size(pop,1)
    x=pop(l,:);
g=0;
for i=2:10
    g=g+(x(i)-(1-0.9*sin(0.2*pi*t)))^2;
end
F1=(1+g)*(x(1)+0.05*sin(w*pi*x(1)));
F2=(1+g)*(s-x(1)+0.05*sin(w*pi*x(1)));
c11=sin(-pi/16)*(F2-1)+cos(-pi/16)*F1;
c1=max(0,0.2*abs(sin(pi*c11)).^0.5-cos(-pi/4)*(F2-1)+sin(-pi/4)*F1+m);
c2=max(0,F1+F2-6);
c3=max(0,m-6);
if c1+c2+c3==0
    feasible_position=[feasible_position;x];
    feasible_cost=[feasible_cost;F1 F2];
else
    infeasible_position=[infeasible_position;x];
end
   end
   feasible_solution=repmat(empty_individual,size(feasible_position,1),1);   
  
   for i=1:size(feasible_position,1)
       feasible_solution(i).position=feasible_position(i,:);
       feasible_solution(i).cost=feasible_cost(i,:)';
   end
    [feasible_solution, F]=NonDominatedSorting(feasible_solution);
    feasible_solution=CalcCrowdingDistance(feasible_solution,F);
    feasible_solution=SortPopulation_2(feasible_solution);
    %%%%%%%%%%%%%%%%搜集等级为1的个体作为引导个体%%%%%%%%%%%%%%%
best_pop=[];
if isempty(F) || isempty(F{1})
    f=pop;
    return;
end
 best_pop1=feasible_solution(F{1});  

 for i=1:size(best_pop1,1)
     best_pop=[best_pop;best_pop1(i).position];
 end
final_pop=[final_pop;best_pop];
 for k=size(best_pop1,1)+1:size(feasible_solution,1)
     infeasible_position=[infeasible_position;[feasible_solution(k).position]] ;
 end
 %%%%%%%%%%%%%%%move%%%%%%%%%%%%%
    for k1=1:size(infeasible_position,1)
     x=infeasible_position(k1,:);
     a=randperm(size(best_pop,1),1);
     x_best=best_pop(a,:);
    new_variable=x+rand*(x_best-x);
       for jj=1:V
               if(new_variable(jj)>max_range(jj))
                    new_variable(jj)=max_range(jj);
               elseif(new_variable(jj)<min_range(jj))
              new_variable(jj)=min_range(jj);
               end
       end
       final_pop=[final_pop;new_variable];
    end
    final_pop_cv=CV(fun,final_pop,t);   
    fea_num= numel(find(final_pop_cv==0));
    f=final_pop;
end
    
    
