function [OP_solution,OP_solution_pf,PF_true,time1,time2,IGD1,IGD2,HV,SP]=main(fun,N,M,V,pps,G,min_range,max_range,options)
if nargin < 9 || isempty(options)
    options=struct();
end
if isfield(options,'seed') && ~isempty(options.seed)
    rng(options.seed);
end

% Default paper settings: N=200, M=2, V=10, pps=100, G=1000.
% The original scripts used the fixed time sequence 0:0.25:20.
if isfield(options,'time_points') && ~isempty(options.time_points)
    m1=options.time_points;
elseif ~isempty(G) && isnumeric(G) && isscalar(G) && G > 0
    m1=0:0.25:min(20,0.25*(G-1));
else
    m1=0:0.25:20;
end
t_end=length(m1);
OP_solution=cell(t_end,1);
OP_solution_pf=cell(t_end,1);
PF_true=cell(t_end,1);
time1=cell(t_end,1);
time2=cell(t_end,1);
IGD1=cell(t_end,1);
IGD2=cell(t_end,1);
HV=cell(t_end,1);
SP=cell(t_end,1);

if isfield(options,'max_batch') && ~isempty(options.max_batch)
    max_batch=options.max_batch;
else
    max_batch=20;  % 批次
end
OP_solution_set=cell(max_batch,1); %用于存放每个批量
OP_solution_pf_set=cell(max_batch,1);
PF_true_set=cell(max_batch,1);        %真实的pf
time1_set=cell(max_batch,1);     % 静态优化使用的时间
time2_set=cell(max_batch,1);      % 从初始化种群到静态优化结束时间
IGD1_set=cell(max_batch,1);%pop对应的IGD
IGD2_set=cell(max_batch,1);%optimalsolution对应的IGD
HV_set=cell(max_batch,1);
SP_set=cell(max_batch,1);


for batch=1:max_batch
    time=0;
for m=m1        % time instance
    variable_list=[];
T=m;
    time=time+1;
    if time<=2     % randomly initialize for second generation
              tic ;%tic2
    t2=clock;
   pop=initialize_variables(fun,N, M, V, min_range,max_range,T);  %初始种群
             pop_pf= initialize_pf(fun,N, M, V,pop,T);
          IGD_measure1= IGD(fun,pop_pf,T,pps);    %% 初始种群的IGD值
          ini_time=toc;
   [F1,last_pop,nsga_time,all_time,F1_last_Gmax]=NSGA2_1(fun,N, M, V,T,pop,pps,time, ini_time,batch,options); % nsga输出结果
        %%%%%%%%%%IGD%%%%%%%%%%%%%%%%%
    PF1=[F1.cost];
            IGD_measure2= IGD(fun,PF1,T,pps) ;          %% ps的IGD值
          HV_measure=hypeIndicatorSampled(PF1');

          if size(PF1,2)==1
              SP_measure=[];
          else
                   SP_measure=spp(PF1',M);
          end
   for i=1:size(F1,1)
       variable_list=[variable_list;F1(i).position];
   end

if time==1
    first_F1=F1;    % t-1时刻
else
    second_F1=F1;   % t时刻
end  


   
    else
        %%%%%%%%%%%%%%%%dynamic response%%%%%%%%%%%%%%
                      tic ;%tic2
    t2=clock;
        pop=dynamic_response_1(fun,N,M,V,first_F1,second_F1,last_pop,T,min_range,max_range,F1_last_Gmax);  % 不是字典
        pop=variable_move(fun,N,M,V,T,min_range,max_range,pop);
                 pop_pf= initialize_pf(fun,N, M, V,pop,T);
          IGD_measure1= IGD(fun,pop_pf,T,pps);    %% 初始种群的IGD值
%    lam=C_NC(fun,N, M, V,T,pop);    %  迭代后的lambda值
%     cost=evaluate_objective(fun,N,M,V,T,pop,lam);   % 转化为无约束问题的目标值
          ini_time=toc;   
  [F1,last_pop,nsga_time,all_time,F1_last_Gmax]=NSGA2_1(fun,N, M, V,T,pop,pps,time, ini_time,batch,options); % nsga输出结果
           %%%%%%%%%%IGD%%%%%%%%%%%%%%%%%
       PF1=[F1.cost];
            IGD_measure2= IGD(fun,PF1,T,pps) ;          %% ps的IGD值
          HV_measure=hypeIndicatorSampled(PF1');

          if size(PF1,2)==1
              SP_measure=[];
          else
                   SP_measure=spp(PF1',M);
          end
          
   for i=1:size(F1,1)
       variable_list=[variable_list;F1(i).position];
   end
          
        first_F1=second_F1;
        second_F1=F1;
        
        
    end
    pof= true_pof(fun,T);             % 真实的pof
    PF_true{time}=pof;
        OP_solution{time}=variable_list;
    OP_solution_pf{time}=PF1;   %保存pareto前沿
        time1{time}=nsga_time;    % 静态优化时间
    time2{time}=all_time;     % 总时间
     IGD1{time}=IGD_measure1;  %% 保存数据,初始种群
       IGD2{time}=IGD_measure2;  %% 保存数据，ps
          HV{time}=HV_measure';  %% 保存数据
             SP{time}=SP_measure;  %% 保存数据

    
end 
PF_true_set{batch}=PF_true;
OP_solution_set{batch}= OP_solution;
OP_solution_pf_set{batch}=OP_solution_pf;
  time1_set{batch}=time1;
  time2_set{batch}=time2;
  IGD1_set{batch}=IGD1;
  IGD2_set{batch}=IGD2;
  HV_set{batch}=HV;
  SP_set{batch}=SP;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%均值%%%%%%%%%%%%%%%%%
OP_solution=OP_solution_set;
OP_solution_pf=OP_solution_pf_set;
PF_true=PF_true_set;
time1=time1_set;
time2=time2_set;
IGD1=IGD1_set;
IGD2=IGD2_set;
HV=HV_set;
SP=SP_set;
end



