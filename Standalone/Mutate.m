%% 变异操作
function y=Mutate(x,mu,sigma,V,xl,xu)


% 多项式变异
x_num=V;
off_1=x;
x_max=xu;
x_min=xl;
yita2=40;
 u2=zeros(1,x_num);
       delta=zeros(1,x_num);
            %r1=rand(1)*fr;
       for j=1:x_num
           u2(j)=rand(1);
           if(u2(j)<0.5)
               delta(j)=(2*u2(j))^(1/(yita2+1))-1;
           else
               delta(j)=1-(2*(1-u2(j)))^(1/(yita2+1));
           end
           off_1(j)=off_1(j)+delta(j);
           %使子代在定义域内
           if(off_1(j)>x_max(j))
               off_1(j)=x_max(j) ;      %+r1*(off_1(j)-x_max(j));
           elseif(off_1(j)<x_min(j))
               off_1(j)=x_min(j);
         
           end
       end
       y=off_1;

end