classdef PWDCMO_DCMOP < PROBLEM
% <2025> <multi> <real> <constrained> <dynamic>
% Dynamic constrained benchmark problems used by PWDCMO
% instance ---  1 --- Problem instance number in 1..8
% taut     --- 20 --- Number of generations for one environment
% nt       ---  4 --- Number of time steps per unit time

%------------------------------- Reference --------------------------------
% D. Gong, M. Rong, N. Hu, Y. Wang, W. Pedrycz, and S. Yang, A prediction
% and weak coevolution-based dynamic constrained multiobjective
% optimization. IEEE Transactions on Evolutionary Computation, 2025,
% 29(4): 1328-1342.
%--------------------------------------------------------------------------

    properties(Access = private)
        instance = 1;
        taut     = 20;
        nt       = 4;
    end
    methods
        function Setting(obj)
            [obj.instance,obj.taut,obj.nt] = obj.ParameterSet(1,20,4);
            obj.instance = min(max(round(obj.instance),1),8);
            obj.taut     = max(1,round(obj.taut));
            obj.nt       = max(1,obj.nt);
            obj.M        = 2;
            if isempty(obj.D); obj.D = 10; end
            obj.lower    = [0,zeros(1,obj.D-1)];
            obj.upper    = [1,2*ones(1,obj.D-1)];
            obj.encoding = ones(1,obj.D);
        end
        function Population = Evaluation(obj,varargin)
            PopDec     = obj.CalDec(varargin{1});
            PopObj     = obj.CalObj(PopDec);
            PopCon     = obj.CalCon(PopDec);
            Population = SOLUTION(PopDec,PopObj,PopCon,zeros(size(PopDec,1),1)+obj.FE);
            obj.FE     = obj.FE + length(Population);
        end
        function PopObj = CalObj(obj,PopDec)
            [PopObj,~] = PWDCMO_DCMOPEvaluate(PopDec,obj.instance,obj.CurrentTime());
        end
        function PopCon = CalCon(obj,PopDec)
            [~,PopCon] = PWDCMO_DCMOPEvaluate(PopDec,obj.instance,obj.CurrentTime());
        end
        function R = GetOptimum(obj,N)
            R = PWDCMO_DCMOPPF(obj.instance,obj.CurrentTime(),N);
        end
        function R = GetPF(obj)
            R = obj.GetOptimum(200);
        end
        function score = CalMetric(obj,metName,Population)
            t      = floor(Population.adds(zeros(length(Population),1))/obj.N/obj.taut)/obj.nt;
            t      = round(t*1e12)/1e12;
            change = [0;find(t(1:end-1)~=t(2:end));length(t)];
            Scores = zeros(1,length(change)-1);
            for i = 1 : length(change)-1
                subPop    = Population(change(i)+1:change(i+1));
                optimum   = PWDCMO_DCMOPPF(obj.instance,t(change(i)+1),10000);
                Scores(i) = feval(metName,subPop,optimum);
            end
            score = mean(Scores);
        end
        function DrawObj(obj,Population)
            t      = floor(Population.adds(zeros(length(Population),1))/obj.N/obj.taut)/obj.nt;
            t      = round(t*1e12)/1e12;
            change = [0;find(t(1:end-1)~=t(2:end));length(t)];
            tempStream = RandStream('mlfg6331_64','Seed',2);
            for i = 1 : length(change)-1
                color = rand(tempStream,1,3);
                Draw(Population(change(i)+1:change(i+1)).objs+(i-1)*0.5,'o','MarkerSize',5,'Markerfacecolor',sqrt(color),'Markeredgecolor',color);
                Draw(PWDCMO_DCMOPPF(obj.instance,t(change(i)+1),200)+(i-1)*0.5,'-','LineWidth',1,'Color',color);
            end
        end
    end
    methods(Access = private)
        function t = CurrentTime(obj)
            t = floor(obj.FE/obj.N/obj.taut)/obj.nt;
        end
    end
end
