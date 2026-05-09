classdef PWDCMO < ALGORITHM
% <2025> <multi> <real> <constrained> <dynamic>
% Prediction and weak coevolution based dynamic constrained MOEA
% keepRatio --- 0.25 --- Ratio of retained solutions after a change
% randRatio --- 0.25 --- Ratio of random solutions after a change

%------------------------------- Reference --------------------------------
% D. Gong, M. Rong, N. Hu, Y. Wang, W. Pedrycz, and S. Yang, A prediction
% and weak coevolution-based dynamic constrained multiobjective
% optimization. IEEE Transactions on Evolutionary Computation, 2025,
% 29(4): 1328-1342.
%------------------------------- Copyright --------------------------------
% This implementation follows the PlatEMO algorithm interface. If this file
% is used inside PlatEMO, publications should also acknowledge PlatEMO.
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            [keepRatio,randRatio] = Algorithm.ParameterSet(0.25,0.25);
            Algorithm.save = sign(Algorithm.save)*inf;
            [W,~] = UniformPoint(Problem.N,Problem.M);

            Population = Problem.Initialization();
            CA = PWDCMO_UpdateCA(Population,Problem.N,W);
            DA = PWDCMO_Select(Population,Problem.N,false);
            [FrontNoCA,CrowdDisCA] = PWDCMO_Fitness(CA,true);
            [FrontNoDA,CrowdDisDA] = PWDCMO_Fitness(DA,false);

            LastPopulation = CA;
            FirstFront     = [];
            SecondFront    = PWDCMO_BestFront(CA);
            LastStaticFront = SecondFront;

            while Algorithm.NotTerminated(CA)
                if PWDCMO_Changed(Problem,CA)
                    OldPopulation = LastPopulation;
                    FirstFront    = SecondFront;
                    CA = Problem.Evaluation(CA.decs);
                    DA = Problem.Evaluation(DA.decs);
                    SecondFront = PWDCMO_BestFront(CA);

                    if ~isempty(FirstFront) && ~isempty(SecondFront)
                        Dec = PWDCMO_DynamicResponse(Problem,OldPopulation,FirstFront,SecondFront,LastStaticFront,keepRatio,randRatio);
                        CA  = Problem.Evaluation(Dec);
                        DA  = CA;
                    else
                        CA = PWDCMO_ChangeFallback(Problem,CA,keepRatio);
                        DA = CA;
                    end
                    LastPopulation = CA;
                    LastStaticFront = PWDCMO_BestFront(CA);
                    CA = PWDCMO_UpdateCA(CA,Problem.N,W);
                    DA = PWDCMO_Select(DA,Problem.N,false);
                    [FrontNoCA,CrowdDisCA] = PWDCMO_Fitness(CA,true);
                    [FrontNoDA,CrowdDisDA] = PWDCMO_Fitness(DA,false);
                end

                nCA = floor(Problem.N/2);
                nCA = nCA + mod(nCA,2);
                nDA = Problem.N - nCA;
                nDA = nDA + mod(nDA,2);
                if nCA + nDA > Problem.N
                    nDA = nDA - 2;
                end

                MatingPoolCA = TournamentSelection(2,nCA,FrontNoCA,-CrowdDisCA);
                MatingPoolDA = TournamentSelection(2,nDA,FrontNoDA,-CrowdDisDA);
                OffspringCA  = OperatorGA(Problem,CA(MatingPoolCA));
                OffspringDA  = OperatorGA(Problem,DA(MatingPoolDA));
                Offspring    = [OffspringCA,OffspringDA];

                CA = PWDCMO_UpdateCA([CA,Offspring],Problem.N,W);
                DA = PWDCMO_Select([DA,Offspring],Problem.N,false);
                [FrontNoCA,CrowdDisCA] = PWDCMO_Fitness(CA,true);
                [FrontNoDA,CrowdDisDA] = PWDCMO_Fitness(DA,false);

                NewFront = PWDCMO_BestFront(CA);
                if ~isempty(NewFront)
                    LastStaticFront = NewFront;
                    SecondFront     = NewFront;
                end
                LastPopulation = CA;
            end
        end
    end
end

function changed = PWDCMO_Changed(Problem,Population)
    if isempty(Population)
        changed = false;
        return;
    end
    sampleNo = max(1,ceil(length(Population)/10));
    RePop1   = Population(randperm(length(Population),sampleNo));
    RePop2   = Problem.Evaluation(RePop1.decs);
    ObjChanged = abs(RePop1.objs-RePop2.objs) > 1e-12;
    ConChanged = abs(RePop1.cons-RePop2.cons) > 1e-12;
    changed  = any(ObjChanged(:)) || any(ConChanged(:));
end

function Dec = PWDCMO_DynamicResponse(Problem,LastPopulation,FirstFront,SecondFront,LastStaticFront,keepRatio,randRatio)
    N     = Problem.N;
    D     = Problem.D;
    nKeep = min(N,floor(N*keepRatio));
    nRand = min(N-nKeep,floor(N*randRatio));
    nPred = N - nKeep - nRand;

    ReLast = Problem.Evaluation(LastPopulation.decs);
    Dec1   = PWDCMO_RetainDec(ReLast,nKeep);
    if nRand > 0
        RandPop = Problem.Initialization(nRand);
        Dec2    = RandPop.decs;
    else
        Dec2 = zeros(0,D);
    end

    FirstDec  = FirstFront.decs;
    SecondDec = SecondFront.decs;
    direction = mean(SecondDec,1) - mean(FirstDec,1);
    if nPred > 0
        baseIndex = mod(0:nPred-1,size(SecondDec,1)) + 1;
        sigma     = PWDCMO_PredictionSigma(SecondFront,LastStaticFront);
        Dec3      = SecondDec(baseIndex,:) + repmat(direction,nPred,1) + sigma.*randn(nPred,D);
    else
        Dec3 = zeros(0,D);
    end

    Dec = Problem.CalDec([Dec1;Dec2;Dec3]);
end

function sigma = PWDCMO_PredictionSigma(SecondFront,LastStaticFront)
    if isempty(LastStaticFront) || isempty(SecondFront)
        sigma = 0;
        return;
    end
    m = min(length(SecondFront),length(LastStaticFront));
    if m == 0
        sigma = 0;
        return;
    end
    F1    = SecondFront(1:m).objs;
    F0    = LastStaticFront(1:m).objs;
    scale = max(max(F1,[],1)-min(F1,[],1),1e-12);
    Delta = abs((F1-F0)./scale);
    sigma = mean(Delta(:));
end

function Population = PWDCMO_ChangeFallback(Problem,Population,keepRatio)
    nKeep = floor(Problem.N*keepRatio);
    Dec1  = PWDCMO_RetainDec(Population,nKeep);
    P2    = Problem.Initialization(Problem.N-size(Dec1,1));
    Population = Problem.Evaluation([Dec1;P2.decs]);
end

function Dec = PWDCMO_RetainDec(Population,N)
    if N <= 0
        Dec = zeros(0,size(Population.decs,2));
        return;
    end
    CV       = sum(max(0,Population.cons),2);
    feasible = CV == 0;
    selected = [];
    if any(feasible)
        FeasiblePop = Population(feasible);
        FeasiblePop = PWDCMO_Select(FeasiblePop,min(N,length(FeasiblePop)),false);
        selected    = FeasiblePop.decs;
    end
    if size(selected,1) < N
        remain = N - size(selected,1);
        infea  = find(~feasible);
        [~,rank] = sort(CV(infea),'ascend');
        infea = infea(rank(1:min(remain,length(rank))));
        selected = [selected;Population(infea).decs];
    end
    Dec = selected;
end

function Front = PWDCMO_BestFront(Population)
    CV       = sum(max(0,Population.cons),2);
    feasible = find(CV == 0);
    if isempty(feasible)
        Front = [];
    else
        FrontNo = NDSort(Population(feasible).objs,1);
        Front   = Population(feasible(FrontNo==1));
    end
end

function Population = PWDCMO_UpdateCA(Population,N,W)
    CV       = sum(max(0,Population.cons),2);
    feasible = CV == 0;
    if sum(feasible) >= N
        Population = PWDCMO_Select(Population(feasible),N,false);
        return;
    end

    CA = Population(feasible);
    remain = N - length(CA);
    Infeasible = Population(~feasible);
    if remain > 0 && ~isempty(Infeasible)
        CVI = CV(~feasible);
        Obj = Infeasible.objs;
        Z   = min(Population.objs,[],1);
        Region = PWDCMO_AssignRegion(Obj,W);
        APD = max(abs(Obj-repmat(Z,size(Obj,1),1))./W(Region,:),[],2);
        Fit = [CVI,APD];
        Next = PWDCMO_SelectByFitness(Fit,remain);
        CA = [CA,Infeasible(Next)];
    end

    if length(CA) < N
        Extra = Population(randperm(length(Population),N-length(CA)));
        CA = [CA,Extra];
    end
    Population = CA;
end

function [FrontNo,CrowdDis] = PWDCMO_Fitness(Population,Constrained)
    if Constrained
        [FrontNo,~] = NDSort(Population.objs,Population.cons,length(Population));
    else
        [FrontNo,~] = NDSort(Population.objs,length(Population));
    end
    CrowdDis = CrowdingDistance(Population.objs,FrontNo);
end

function Population = PWDCMO_Select(Population,N,Constrained)
    if isempty(Population) || N <= 0
        Population = Population([]);
        return;
    end
    N = min(N,length(Population));
    if Constrained
        [FrontNo,MaxFNo] = NDSort(Population.objs,Population.cons,N);
    else
        [FrontNo,MaxFNo] = NDSort(Population.objs,N);
    end
    CrowdDis = CrowdingDistance(Population.objs,FrontNo);
    Next     = FrontNo < MaxFNo;
    Last     = find(FrontNo == MaxFNo);
    [~,Rank] = sort(CrowdDis(Last),'descend');
    Next(Last(Rank(1:N-sum(Next)))) = true;
    Population = Population(Next);
end

function Next = PWDCMO_SelectByFitness(Fit,N)
    N = min(N,size(Fit,1));
    [FrontNo,MaxFNo] = NDSort(Fit,N);
    CrowdDis = CrowdingDistance(Fit,FrontNo);
    Next = FrontNo < MaxFNo;
    Last = find(FrontNo == MaxFNo);
    [~,Rank] = sort(CrowdDis(Last),'descend');
    Next(Last(Rank(1:N-sum(Next)))) = true;
end

function Region = PWDCMO_AssignRegion(PopObj,W)
    PopNorm = sqrt(sum(PopObj.^2,2));
    WNorm   = sqrt(sum(W.^2,2))';
    Cosine  = PopObj*W'./max(PopNorm*WNorm,eps);
    [~,Region] = max(Cosine,[],2);
end
