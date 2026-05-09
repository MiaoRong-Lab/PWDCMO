function R = PWDCMO_DCMOPPF(instance,t,N)
% Approximate the true Pareto front of a PWDCMO DCMOP instance.

    if nargin < 3
        N = 10000;
    end
    [s,m,w] = localParameters(instance,t);
    sampleN = max(N,2001);

    x1 = linspace(0,1,sampleN)';
    f1 = x1 + 0.05*sin(w*pi*x1);
    f2 = s - x1 + 0.05*sin(w*pi*x1);
    optDec = [x1,repmat(1-0.9*sin(0.2*pi*t),sampleN,9)];
    [~,curveCon] = PWDCMO_DCMOPEvaluate(optDec,instance,t);
    feasibleCurve = all(curveCon <= 0,2);
    P = [f1(feasibleCurve),f2(feasibleCurve)];

    a = 0:0.005:4;
    b = 0:0.005:6;
    [A,B] = meshgrid(a,b);
    c11 = sin(-pi/16).*(B-1) + cos(-pi/16).*A;
    c1  = 0.2*abs(sin(pi*c11)).^0.5 - cos(-pi/4).*(B-1) + sin(-pi/4).*A + m;
    feasible = c1 <= 0 & A+B-6 <= 0 & m-6 <= 0;
    hasFeasible = any(feasible,1);
    if any(hasFeasible)
        bMin = nan(1,length(a));
        for i = find(hasFeasible)
            bColumn = B(:,i);
            bMin(i) = min(bColumn(feasible(:,i)));
        end
        [curveF1,order] = sort(f1);
        curveF2 = f2(order);
        [curveF1,uniqueIndex] = unique(curveF1,'stable');
        curveF2 = curveF2(uniqueIndex);
        query = a(hasFeasible);
        lowerBoundary = bMin(hasFeasible);
        curveAtQuery = interp1(curveF1,curveF2,query,'linear','extrap');
        boundary = [query(:),max(curveAtQuery(:),lowerBoundary(:))];
        boundary = boundary(all(isfinite(boundary),2),:);
        P = [P;boundary];
    end

    if isempty(P)
        R = [f1,f2];
    else
        R = localNondominated(P);
    end
    if size(R,1) > N
        index = unique(round(linspace(1,size(R,1),N)));
        R = R(index,:);
    end
end

function P = localNondominated(P)
    P = unique(round(P*1e12)/1e12,'rows');
    P = sortrows(P,[1,2]);
    keep = false(size(P,1),1);
    best = inf;
    for i = 1 : size(P,1)
        if P(i,2) < best - 1e-12
            keep(i) = true;
            best = P(i,2);
        end
    end
    P = P(keep,:);
end

function [s,m,w] = localParameters(instance,t)
    switch instance
        case 1
            s = max(3.5-0.14*t,0.7+0.14*t);
            m = max(1.43-0.05*t,0.43+0.05*t);
            w = 2;
        case 2
            s = max(2.5-0.05*t,1.5+0.05*t);
            m = max(1.16-0.075*t,-0.43+0.075*t);
            w = 2;
        case 3
            s = min(2.1+0.14*t,4.9-0.14*t);
            m = min(0.93+0.05*t,1.93-0.05*t);
            w = 2;
        case 4
            s = min(2+0.05*t,3-0.05*t);
            m = min(0.41+0.075*t,1.91-0.075*t);
            w = 2;
        case 5
            s = max(3.5-0.14*t,0.7+0.14*t);
            m = max(1.43-0.05*t,0.43+0.05*t);
            w = 6*sin(0.2*pi*(t+1));
        case 6
            s = max(2.5-0.05*t,1.5+0.05*t);
            m = max(1.16-0.075*t,-0.43+0.075*t);
            w = 6*sin(0.2*pi*(t+1));
        case 7
            s = min(2.1+0.14*t,4.9-0.14*t);
            m = min(0.93+0.05*t,1.93-0.05*t);
            w = 6*sin(0.2*pi*(t+1));
        otherwise
            s = min(2+0.05*t,3-0.05*t);
            m = min(0.41+0.075*t,1.91-0.075*t);
            w = 6*sin(0.2*pi*(t+1));
    end
end
