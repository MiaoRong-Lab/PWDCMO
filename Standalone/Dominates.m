%% ÅĞ¶ÏÖ§Åä¹ØÏµ
function b=Dominates(x,y)
    if isstruct(x)
        x=x.cost;
    end
    if isstruct(y)
        y=y.cost;
    end
    b=all(x<=y) & any(x<y);
end
