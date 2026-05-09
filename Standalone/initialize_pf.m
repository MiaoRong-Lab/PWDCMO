function f = initialize_pf(fun,N, M, V,pop,t)

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


pf=[];
for i=1:size(pop,1)
    x=pop(i,:);
    g=0;
for i1=2:10
    g=g+(x(i1)-(1-0.9*sin(0.2*pi*t)))^2;
end
f1=(1+g)*(x(1)+0.05*sin(w*pi*x(1)));
f2=(1+g)*(s-x(1)+0.05*sin(w*pi*x(1)));
pf=[pf [f1;f2]];
end
f=pf;
end