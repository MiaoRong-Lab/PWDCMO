function f = true_pof(fun,t)



%%%%%%%%%%%%%%%%%%论文中的pof%%%%%%%%%%%%%%%%%%%%%%
  
sample=[];
pof=[];
pof1=[];
noncontraint_pop=[];
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

% 生成pos
% x1=linspace(0,1,500);
x1=(0:0.0005:1);
% for i=2:10
%     x(i)=1-0.9*sin(0.2*pi*t);
%     g=g+(x(i)-(1-0.9*sin(0.2*pi*t)))^2;
% end
% pof计算
F1=(x1+0.05*sin(w*pi*x1));
F2=(s-x1+0.05*sin(w*pi*x1));



     
      %%%%%%%%%%%%可行域%%%%%%%%%%%%%%%%%%%%%%5
     sample2=[];
      a1=0;
      a=(0:0.005:4);      % 离散范围
      b=(0:0.005:6);
      for k=1:numel(a)        % 均匀采点
          y11=a(k);
          for i=1:numel(b)
              y22=b(i);
              c11=sin(-pi/16)*(y22-1)+cos(-pi/16)*y11;
               c1=max(0,0.2*abs(sin(pi*c11))^0.5-cos(-pi/4)*(y22-1)+sin(-pi/4)*y11+m);
               % c1=max(0,0.2*abs(sin(pi*c11)).^0.5-cos(-pi/4)*(f2-1)+sin(-pi/4)*f1+m);
               c2=max(0,y11+y22-6);
                c3=max(0,m-6);
             if c1+c2+c3==0
              sample=[sample [y11;y22]];          % 在区域中采点
               end
          end
      end

      sample1=sample(1,:);
      k=find(sample1>F1(numel(F1)));
      if ~isempty(k)
          sample(:,(k(1):size(sample,2)))=[];
      end
      

       
      for ii=1:sample(1,size(sample,2))/0.005+1
          sample1=sample(1,:);
          a1_index=find(abs(sample1-a1)<1e-5); 
            a2_index=find(abs(F1-a1)<0.0005);
            if isempty(a2_index)
                break;
            end
            a2_index=a2_index(1);
          noncontraint_pop=[noncontraint_pop [F1(a2_index);F2(a2_index)]];  
              F1(:,1:a2_index)=[];
              F2(:,1:a2_index)=[];
         k=length(a1_index); 
         sample2=[sample2 sample(:,1)];
         sample(:,(1:k))=[];
         a1=a1+0.005;
      end
      
          for d=1:size( noncontraint_pop,2)
  f1=noncontraint_pop(1,d);
  f2=noncontraint_pop(2,d);
c11=sin(-pi/16)*(f2-1)+cos(-pi/16)*f1;
c1=max(0,0.2*abs(sin(4*pi*c11)).^0.5-cos(-pi/4)*(f2-1)+sin(-pi/4)*f1+m);
 c2=max(0,f1+f2-6);
c3=max(0,m-6);
if c1+c2+c3==0
pof=[pof [f1;f2]];                        %% 满足约束的无约束pof
end
          end
      


      
      %%%%%%%%%%%非支配排序%%%%%%%%%%%%%%%
  for i3=1:size(sample2,2)
         value1=sample2(2,i3);
         value2=noncontraint_pop(2,i3);
         if value1>value2
             pof1=[pof1 sample2(:,i3)];
         else
             pof1=[pof1 noncontraint_pop(:,i3)];
         end
     end
       empty.cost=[];    
 empty.Rank=[];
 pop=repmat(empty,size(pof1,2),1);   % 复制
  for i22=1:size(pof1,2)
    pop(i22).cost=pof1(:,i22);
 end
  [pop, F]=NonDominatedSorting(pop);        %% 非支配排序
  F1=pop(F{1});
  POF=[F1.cost];



      f=POF;
end
      




