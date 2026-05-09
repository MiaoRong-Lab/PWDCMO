function [ff1,ff2,ff3,ff4,ff5] = NSGA2_1(fun,NN, Mm, VV,t,popp,pps,time,ini_time,batch,options)
if nargin < 11 || isempty(options)
    options=struct();
end
if ~isfield(options,'plot')
    options.plot=false;
end
if ~isfield(options,'save_figures')
    options.save_figures=false;
end
if ~isfield(options,'figure_dir')
    options.figure_dir=fullfile(pwd,'results','figures');
end
if ~isfield(options,'verbose')
    options.verbose=false;
end
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

tic;
%% 初始化
T=t;
N=NN;
M=Mm;
V=VV;
VarMin=[0 0 0 0 0 0 0 0 0 0];          % 下限
VarMax= [1 2 2 2 2 2 2 2 2 2];          % 上限
step=[1,1,1,1,1,1,1,1,1,1];
len=(VarMax-VarMin)./step;
var=[VarMin;step;VarMax;round(len,0)];
if time<=2
Gmax=80;           %迭代次数
else
    Gmax=50;
end
nPop=N;    % 种群个体个数
pCrossover=0.8;                         % 交叉率
nCrossover=2*round(pCrossover*nPop/2);  % 子代个数
pMutation=0.1;                          % 变异率
nMutation=round(pMutation*nPop);        % 变异数量
mu=0.02;                    % 突变概率
sigma=0.1*(VarMax-VarMin);  % 突变步长
%%%%%weight vector
 [W,N] = UniformPoint(size(popp,1),M);
 CA=popp;
 DA=popp;

%% 问题定义
empty_individual.position=[];     
empty_individual.cost=[];        
empty_individual.Rank=[];
empty_individual.DominationSet=[];
empty_individual.DominatedCount=[];
empty_individual.Fitness=[];

 %%%%%%%%%%%%%%%%%%%DA不考虑约束违反%%%%%%%%%%%%%%%%
 DA_pf=initialize_pf(fun,N, M, V,DA,t);
 pop_DA=repmat(empty_individual,size(DA,1),1);   % 复制
 for i=1:size(DA,1)
    pop_DA(i).position=DA(i,:);
    pop_DA(i).cost=DA_pf(:,i);
 end
    % 非支配排序
[pop_DA, F]=NonDominatedSorting(pop_DA);
% 计算拥挤度
pop_DA=CalcCrowdingDistance(pop_DA,F);   % fitness里面放拥挤距离
% 对个体排序
[pop_DA, F]=SortPopulation_2(pop_DA);
pop_nc=pop_DA(1:min(pps,numel(pop_DA)));
F1=pop_nc;
F1_last_Gmax=pop_nc;

%% NSGA-II 主循环
for it=1:Gmax
% 计算约束违反CA
CA_pf=initialize_pf(fun,N, M, V,CA,t);
CA_cv=CV(fun,CA,t);
% %% % 找到可行解的索引
 current_ca=find(CA_cv==0);
 pop_ca=repmat(empty_individual,size(current_ca,1),1);   % 复制
 for i=1:length(current_ca)
     pop_ca(i).position=CA(current_ca(i),:);
     pop_ca(i).cost=CA_pf(:,current_ca(i));
 end
     % CA中可行解非支配排序
[pop_ca, F]=NonDominatedSorting(pop_ca);
% 计算拥挤度
pop_ca=CalcCrowdingDistance(pop_ca,F);   
% 对个体排序
[pop_ca, F]=SortPopulation_2(pop_ca);
last_rank1=pop_ca(length(current_ca)).Rank;
if length(current_ca)<length(CA_cv)

  chromosome_infea=[];
current_c1=find(CA_cv>0);   
for i=1:length(current_c1)
    chromosome_infea=[chromosome_infea;CA(current_c1(i),:) [CA_pf(:,current_c1(i))]' last_rank1+1 CA_cv(current_c1(i))];
end

chromosome_infea=sortrows(chromosome_infea,size(chromosome_infea,2));

for i=1:size(chromosome_infea,1)
    pop_ca(length(current_ca)+i).position=chromosome_infea(i,1:V);
    pop_ca(length(current_ca)+i).cost=chromosome_infea(i,V+1:V+M)';
    pop_ca(length(current_ca)+i).Rank=last_rank1+1;
    pop_ca(length(current_ca)+i).Fitness=chromosome_infea(i,end);
end
end
 
Q=[];   
    for k=1:nCrossover/2     
     p1=tournamentsel(pop_ca,last_rank1+1);
     p2=tournamentsel(pop_ca,last_rank1+1);
     while isequal(p2.position,p1.position)
         p2=tournamentsel(pop_ca,last_rank1+1);
     end   
    [offspring1, offspring2]=Crossover(p1.position,p2.position,V,VarMin,VarMax) ;
        Q=[Q;offspring1;offspring2];
    end
        
    for kk=1:nMutation/2
        i=randi([1 size(CA,1)]) ;   
        p=CA(i,:);          
        offspring=Mutate(p,mu,var,V,VarMin,VarMax);
        Q=[Q;offspring];
    end

    for k=1:nCrossover/2     
     p1=tournamentsel_1(pop_DA);
     p2=tournamentsel_1(pop_DA);
     while isequal(p2.position,p1.position)
         p2=tournamentsel_1(pop_DA);
     end   
    [offspring1, offspring2]=Crossover(p1.position,p2.position,V,VarMin,VarMax) ;
        Q=[Q;offspring1;offspring2];
    end
    % 变异
    for kk=1:nMutation/2
       i=randi([1 size(pop_DA,1)]) ;   
        p=pop_DA(i,:);
        offspring=Mutate(p.position,mu,var,V,VarMin,VarMax);
        Q=[Q;offspring];
    end
CA1=[CA CA_pf' CA_cv];
Q_cv=CV(fun,Q,t);
Q_pf=initialize_pf(fun,N, M, V,Q,t);
Q1=[Q Q_pf' Q_cv];

CA = UpdateCA(CA1,Q1,V,M,W);

for i=1:size(Q,1)
    pop_DA(N+i).position=Q(i,:);
     pop_DA(N+i).cost=Q_pf(:,i);
end
 % 非支配排序
[pop_DA, F]=NonDominatedSorting(pop_DA);
% 计算拥挤度
pop_DA=CalcCrowdingDistance(pop_DA,F);   % fitness里面放拥挤距离
% 对个体排序
[pop_DA, F]=SortPopulation_2(pop_DA);
% 防止溢出
pop_DA=pop_DA(1:N);

        %%% 计算当前种群的可行非支配解
        population3=[];
        nc_obc=[];
     for k1=1:size(CA,1)
            individual= CA(k1,:);
             g=0;
              for ii=2:10
                      g=g+( individual(ii)-(1-0.9*sin(0.2*pi*T)))^2;
              end
                  f1=(1+g)*( individual(1)+0.05*sin(w*pi* individual(1)));
                  f2=(1+g)*(s- individual(1)+0.05*sin(w*pi* individual(1)));
                  c11=sin(-pi/16)*(f2-1)+cos(-pi/16)*f1;
                  c1=max(0,0.2*abs(sin(pi*c11))^0.5-cos(-pi/4)*(f2-1)+sin(-pi/4)*f1+m);
                  c2=max(0,f1+f2-6);
                  c3=max(0,m-6);
                        if c1+c2+c3==0
                        
                          population3=[ population3;individual];      
                             nc_obc=[nc_obc [f1;f2]];                 
                        end

     end
     if size(population3,1)>0
             pop_nc=repmat(empty_individual,size(population3,1),1);
   for i1=1:size(population3,1)
   pop_nc(i1).position=population3(i1,:);
   pop_nc(i1).cost= nc_obc(:,i1);
   end
   % 非支配排序
[ pop_nc, F]=NonDominatedSorting( pop_nc);
% 计算拥挤度
pop_nc=CalcCrowdingDistance(pop_nc,F);
% 对个体排序
[pop_nc, F]=SortPopulation_2(pop_nc);
         F1=pop_nc(F{1}); 
       

%     %结果显示
    if options.verbose
        disp(['Iteration ' num2str(it) ': Number of F1 Members = ' num2str(numel(F1))]);
    end
     end  
     if it==Gmax-1
         F1_last_Gmax=pop_nc(1:min(pps,numel(pop_nc)));
     end


     end  
nsga_time=toc;%% 静态优化时间
all_time=nsga_time+ini_time;
%save time nsga_time;
               if options.verbose
                   disp(['静态优化时间：',num2str(nsga_time),'总时间：',num2str(all_time)]);
               end
%     % 动态画图

pof= true_pof(fun,T);             % 真实的pof
if options.plot
    if size(F1,1)>pps
        F2=F1(1:pps);
    else
        F2=F1;
    end
    cost=[F2.cost];       % 绘制PS
    y1=cost(1,:);
    y2=cost(2,:);
    f1=pof(1,:);
    f2=pof(2,:);

    if time==1
        figure;
    end
    sz=5;
    scatter(f1+0.5*T,f2+0.5*T,sz,'r')
    hold on
    sz=10;
    scatter(y1+0.5*T,y2+0.5*T,sz,'*b')
    hold on
    xlabel('f1+0.5*t')
    ylabel('f2+0.5*t')
    legend('PF','DCMOEA');
    title('Ins1')
    drawnow;
    if options.save_figures && time==81
        if ~exist(options.figure_dir,'dir')
            mkdir(options.figure_dir);
        end
        saveas(gcf,fullfile(options.figure_dir,[num2str(batch),'.fig']));
    end
end
   ff1=pop_nc(1:min(pps,numel(pop_nc)));
    ff2=CA;
         ff3=nsga_time;
     ff4=all_time;
     ff5=F1_last_Gmax;

end
    
    


    
        
        
        
        
