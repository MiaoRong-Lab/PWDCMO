
function f = IGD(fun,pf,t,pps)
true_pf= true_pof(fun,t);
pf_number=size(true_pf,2) ;      
sum_distance=0;
for i=1:pf_number
    y1=true_pf(:,i);           
    distance=[];
    for i1=1:size(pf,2)
        y2=pf(:,i1);
        dis=sqrt((y1(1)-y2(1))^2+(y1(2)-y2(2))^2);
        distance=[distance dis];
    end
    sum_distance= sum_distance+min(distance);
end
IGD=sum_distance/pf_number;
f=IGD;
end
    
