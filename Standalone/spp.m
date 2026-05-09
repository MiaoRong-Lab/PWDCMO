function f=spp(X,M)  

D=[];
 for i=1:1:size(X,1)    
     d=0;dd=[];
     y1=X(i,:);  
     for j=1:1:size(X,1)
         if(j~=i)    
             y2=X(j,:);
             ddd=0;
             for jj=1:M  
                 ddd=ddd+(y1(jj)-y2(jj))^2;   
               %ddd=ddd+abs(X(j,k+jj)-X(i,k+jj));   
             end
            dd=[dd;sqrt(ddd)];
         end
     end
    d=min(dd);  
    D=[D;d];    
end
avd=mean(D);%计算平均值
sumd=0;
for i=1:1:size(X,1)
    sumd=sumd+(avd-D(i))^2;
end
f=sqrt(sumd/(size(X,1)-1));    
end
    