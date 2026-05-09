%% ½»²æ²Ù×÷
function [y1, y2]=Crossover(x1,x2,V,xl,xu)


        
        % SBX½»²æ
        x_num=V;
        chromo_parent_1=x1;
        chromo_parent_2=x2;
        x_max=xu;
        x_min=xl;
        yita1=5;
        u1=zeros(1,x_num);
       gama=zeros(1,x_num);

       for j=1:x_num
           u1(j)=rand(1);
           if u1(j)<0.5
               gama(j)=(2*u1(j))^(1/(yita1+1));
           else
               gama(j)=(1/(2*(1-u1(j))))^(1/(yita1+1));
           end
           off_1(j)=0.5*((1+gama(j))*chromo_parent_1(j)+(1-gama(j))*chromo_parent_2(j));
           off_2(j)=0.5*((1-gama(j))*chromo_parent_1(j)+(1+gama(j))*chromo_parent_2(j));
         
           if(off_1(j)>x_max(j))
               off_1(j)=x_max(j) ;     %+r*(off_1(j)-x_max(j));
           elseif(off_1(j)<x_min(j))
               off_1(j)=x_min(j);
               %+r*(off_1(j)-x_min(j));
           end
           if(off_2(j)>x_max(j))
               off_2(j)=x_max(j);    %+r*(off_2(j)-x_max(j));
           elseif(off_2(j)<x_min(j))
               off_2(j)=x_min(j);
               %+r*(off_2(j)-x_min(j));
           end
       end
       y1=off_1;
       y2=off_2;

        


end