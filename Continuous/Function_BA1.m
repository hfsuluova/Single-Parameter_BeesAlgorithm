%___________________________________________________________________%
%  The Single-parameter Bees Algorithm (BA1) source codes           %
%  version 1.0                                                      %
%                                                                   %
%  Developed in MATLAB R2022a                                       %
%                                                                   %
%  Programmer: Hamid Furkan Suluova                                 %
%  Authors: H.F. Suluova and D.T. Pham                              %
%         e-Mail: hfsuluova@gmail.com                               %
%                 hfs158@bham.ac.uk                                 %
%___________________________________________________________________%

function [it, OptCost, NFE] = Function_BA1(F, n, MaxEval)
%% Problem Definition

[lb,ub,dim,fobj] = Function_Library(F);

Dims = dim;
ObjFunction = fobj; % Objective Function
VarSize = [1 Dims]; % Decision Variables Matrix Size
VarMin = lb; % Decision Variables Lower Bound
VarMax = ub; % Decision Variables Upper Bound
range = VarMax-VarMin;

%% Parameters
%n=100;
%MaxEval = 500000; 
MaxIt = round(MaxEval/n);
%MaxIt = 1000;

%% Initialization

Unknown_Bee.Name = [];
Unknown_Bee.Position = [];
Unknown_Bee.Cost = [];
Unknown_Bee.Size = [];
Unknown_Bee.Stagnated = [];
Unknown_Bee.counter = [];
Unknown_Bee.Cluster = [];
Bee = repmat(Unknown_Bee,n,1);
counter = 0;

% Distributing Bees to Search Space
for i = 1:n
    Bee(i).Name = i;
    Bee(i).Position = unifrnd(VarMin,VarMax,VarSize);
    Bee(i).Cost = ObjFunction(Bee(i).Position);
    Bee(i).Size = range;
    Bee(i).Stagnated = 0;
    Bee(i).Cluster = 0;
    counter = counter+1;
    Bee(i).counter = counter;
end


%% Bee Sorting 
[~, RankOrder] = sort([Bee.Cost]);
Bee = Bee(RankOrder);

%% Euclidean Distance to Best Bee
BestBee = Bee(1);
for i=1:n
    sq_sum = 0;
    for j = 1:Dims
        sq_sum = sq_sum + (BestBee.Position(j) - Bee(i).Position(j))^2;
    end
    Bee(i).Distance = sqrt(sq_sum);
end

 
%% Incremental K-means

B=[Bee.Name; Bee.Distance].';
for K = 1:(n-1)
    [idx,C,~,D] = kmeans(B(:,2),K);
    Clustering(K).Name=K;
    Clustering(K).Indexes= idx;
    Clustering(K).Centroids = C;
    Clustering(K).Distances = D;
    Dist_Sq_Sum = zeros(K,1);

    for i=1:K
        for j=1:n
            if ismember(D(j,i),D(idx==i,i))
                Dist_Sq_Sum(i) = Dist_Sq_Sum(i) + (C(i)-D(j,i))^2;
                Clustering(K).Distortion(i) = Dist_Sq_Sum(i);
            end
        end
    end
    Clustering(K).Sum_Distortion = sum(Dist_Sq_Sum);
    
%Calculation alpha_K and S_k
    if (K==2) && Dims>1
        alpha(K) = 1 - (3/(4*Dims));
    elseif K>2 && Dims>1
        alpha(K) = alpha(K-1) + ((1 - alpha(K-1))/6);
    else
        alpha(K) = 1;
    end
    
    if K>1 && Clustering(K-1).Sum_Distortion~=0
        Clustering(K).f_K= Clustering(K).Sum_Distortion/(alpha(K)*Clustering(K-1).Sum_Distortion);
    else
        Clustering(K).f_K=1;
    end
end

for l=1:n-1
    Res(l)=Clustering(l).Sum_Distortion;
    Res2(l)=Clustering(l).f_K;
end
[min_dist,min_dist_idx] = min(Res,[],'linear');
[SelFunc,SelFuncidx] = min(Res2,[],'linear');


% min_dist_idx is the output based on the Sum of the cluster distortion
% SelFuncidx is the output based on the f(K) selection function

% min_dist_idx = SelFuncidx; % Can change if desire


R(1)=1;
for x = 2:n-1
    R(x) = Clustering(x).Sum_Distortion/Clustering(x-1).Sum_Distortion;
end



for i=1:n
    Bee(i).Cluster = Clustering(min_dist_idx).Indexes(i);
end

%% Recruitment

% min_dist_idx is the output of clustering process
% n is the number of bees in the colony

% Sorting according to Cluster


for nr=1:min_dist_idx
    Bee_Recruitment(nr) = sum(Clustering(min_dist_idx).Indexes == nr);
end

[~, RecSort]=sort(Bee_Recruitment,'descend');
Bee_Recruitment = Bee_Recruitment(RecSort);


for i = 1:min_dist_idx
    Cluster_Bee = 1;
    for j = 1:n
        if Bee(j).Cluster == i
            Patches(i).Bee(Cluster_Bee) = Bee(j);
            Cluster_Bee = Cluster_Bee + 1; 
        end
    
    end
    Patch(i).Position = Patches(i).Bee(1).Position;
    Patch(i).Cost = Patches(i).Bee(1).Cost;
    Patch(i).Size = Patches(i).Bee(1).Size;
    Patch(i).Stagnated = Patches(i).Bee(1).Stagnated;
    Patch(i).counter = Patches(i).Bee(1).counter;
    Patch(i).Recruited = Bee_Recruitment(i);

end


%% Incremental K Means Figures
% 
% figure(1); plot(R,'-o'); title('Sk/Sk-1'); axis([0 inf 0 2]); grid on
% figure(2); plot(Res); title('Distortion'); grid on
% figure(3); plot(Res2); title('f_K'); axis([0 inf 0 2]); grid on
% figure(4); plot(alpha); title('alpha'); axis([0 inf 0 2]); grid on
% 

%% Bees Algorithm Local and Global Search
BestSol.Cost = inf;
 ssize = linspace(0,1,min_dist_idx);
for it = 1:MaxIt
    if counter >= MaxEval
        break;
    end
    for i = 1:min_dist_idx
        bestWorker.Cost = inf;
        assigntment = D_Triangular (0,ssize(i),1,1,Bee_Recruitment(i));
        for j = 1:Bee_Recruitment(i)
            Worker.Position = Foraging (Patch(i).Position,assigntment(j),VarMax,VarMin,Patch(i).Size);
            Worker.Cost = ObjFunction(Worker.Position);
            Worker.Size = Patch(i).Size;
            Worker.Stagnated = Patch(i).Stagnated;
            Worker.Recruited = Patch(i).Recruited;
            counter = counter+1;
            Worker.counter = counter;
            if Worker.Cost < bestWorker.Cost
                bestWorker = Worker;
            end
        end
        if bestWorker.Cost < Patch(i).Cost
            Patch(i) = bestWorker;
            Patch(i).Stagnated = 0;
        else
            Patch(i).Stagnated = Patch(i).Stagnated+1;
            Patch(i).Size = Patch(i).Size*(1-(3*it/(4*MaxIt)));
            if  Patch(i).Stagnated>=round(n/min_dist_idx)
                Patch(i).Position = unifrnd(VarMin,VarMax,VarSize);
                Patch(i).Cost = ObjFunction(Patch(i).Position);
                Patch(i).Size = range;
                Patch(i).Stagnated = 0;
                counter = counter+1;
                Patch(i).counter = counter;
            end
        end
    end
    % SORTING
    [~, RankOrder] = sort([Patch.Cost]);
    Patch = Patch(RankOrder);
    % Update Best Solution Ever Found
    OptSol = Patch(1);
    if OptSol.Cost < BestSol.Cost
        BestSol = OptSol;
    end
    % taking of result
    OptCost(it) = BestSol.Cost;
    Counter(it) = counter;
    %Time(it) = toc;
    NFE = counter;
    % Display Iteration Information
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(OptCost(it)) '; Fittness Evaluations = ' num2str(Counter(it))]);
    
end
%% Results
% figure;
% semilogy(OptCost,'LineWidth',2);
% xlabel('Iteration');
% ylabel('Best Cost');
end
