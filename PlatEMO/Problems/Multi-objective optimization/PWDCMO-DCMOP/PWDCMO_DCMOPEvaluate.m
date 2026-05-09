function [PopObj,PopCon] = PWDCMO_DCMOPEvaluate(PopDec,instance,t)
% Objective and constraint functions of the PWDCMO DCMOP instances.

    [s,m,w] = PWDCMO_DCMOPParameters(instance,t);
    g       = sum((PopDec(:,2:end) - (1-0.9*sin(0.2*pi*t))).^2,2);
    wave    = 0.05*sin(w*pi*PopDec(:,1));
    f1      = (1+g).*(PopDec(:,1) + wave);
    f2      = (1+g).*(s - PopDec(:,1) + wave);
    PopObj  = [f1,f2];

    c11 = sin(-pi/16).*(f2-1) + cos(-pi/16).*f1;
    c1  = 0.2*abs(sin(pi*c11)).^0.5 - cos(-pi/4).*(f2-1) + sin(-pi/4).*f1 + m;
    c2  = f1 + f2 - 6;
    c3  = zeros(size(f1)) + m - 6;
    PopCon = [c1,c2,c3];
end

function [s,m,w] = PWDCMO_DCMOPParameters(instance,t)
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
