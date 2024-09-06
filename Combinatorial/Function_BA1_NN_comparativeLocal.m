function [sonuc, counter, OptCost] = Function_BA1_NN_comparativeLocal (typeOfFunction,n)
%% Problem Definition
%[typeOfFunction] = 'Eil51'; %{'A280','Att532','Berlin52','Eil51','Eil76',
% 'Fl1577','KroA100','KroA150','KroA200','KroB100','KroB150','KroB200',
% 'KroC100','KroD100','KroE100','Lin318','Pcb442','Pr76','Rat99','Rat783',
% 'St70',Fl417, Tsp225, D198, Ch150,Pr107}
Instance = Tsplib(typeOfFunction);
Dims = Instance.dim;
ObjFunction = @(x) Instance.evaluation(x);% Objective Function
VarSize = [1 Dims];                       % Decision Variables Matrix Size
NNMDist = Instance.D;

%% Bees Algorithm Parameters
MaxEval = 10000000;
%n = 800;
assigntment = linspace(0,1,n);
MaxIt = 1000;    % Maximum Number of Iterations
MaxIt = round(MaxEval/n);
%% Initialization
Unknown_Patch.Position = [];
Unknown_Patch.Cost = [];
Unknown_Patch.Cluster = [];
Unknown_Patch.counter = [];

Scout = repmat(Unknown_Patch,n,1);
counter = 0;
%% Generate Initial Solutions
for i = 1:n
    Scout(i).Position = randperm(Dims);
    Scout(i).Cost = ObjFunction(Scout(i).Position);
    counter = counter+1;
    Scout(i).counter = counter;
end

%% Nearest Neighbourhood
% 
% for i = 1:n
%     s = randi(Dims);
%     IS=[];
%     IS= [IS, s];
%     DistList_Changing = NNMDist;
%     for j = 1:Dims
%         s = IS(end);
%       DistList_Changing(s,:) = inf;
%         temp= DistList_Changing(:,s);
%         [~, minIDX] = min(temp);
%         IS = [IS minIDX];
%     end
%     Scout(i).Position = IS;
%     Scout(i).Cost = ObjFunction(Scout(i).Position);
%     counter = counter+1;
%     Scout(i).counter = counter;
% end
% size = linspace(1,1,n);

%% Sites Selection 
[~, RankOrder] = sort([Scout.Cost]);
Scout = Scout(RankOrder);

%% Euclidean Distance to Best Bee
BestBee = Scout(1);

for i=1:n
    sq_sum = 0;
    for j = 1:Dims
        sq_sum = sq_sum + (BestBee.Position(j) - Scout(i).Position(j));
    end
    Scout(i).Distance = sq_sum;
end

%% Preparing Position to clustering

 B=[Scout.Distance].';


%% Incremental K-means
if n > 50
    lmt = 50;
else 
    lmt = n;
end
for K = 1:lmt
    %[idx,C,~,D] = kmeans(A(:,2:end),K);
    %[idx1,C1,~,D1] = kmeans(A1(:,2:(Dims+1)),K);
    [idx,C,~,D] = kmeans(B,K);
    %[idx,C,~,D] = kmeans(X(:,2),K);
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


for l=1:lmt
    Res(l)=Clustering(l).Sum_Distortion;
    Res2(l)=Clustering(l).f_K;
end
[min_dist,min_dist_idx] = min(Res,[],'linear');
[val,validx] = min(Res2,[],'linear');


for i=1:n
    Scout(i).Cluster = Clustering(min_dist_idx).Indexes(i);
end

%% Recruitment
% min_dist_idx is the output of clustering process
% n is the number of bees in the colony

%Sorting according to Cluster

[~, RankOrder]=sort([Scout.Cost]);
Scout = Scout(RankOrder);


for nr=1:min_dist_idx
    Bee_Recruitment(nr) = sum(Clustering(min_dist_idx).Indexes == nr);
end

for i = 1:min_dist_idx
    Cluster_Bee = 1;
    for j = 1:n
        if Scout(j).Cluster == i
            Patches(i).Scout(Cluster_Bee) = Scout(j);
            Cluster_Bee = Cluster_Bee + 1; 
        end
    
    end
    Patch(i).Position = Patches(i).Scout(1).Position;
    Patch(i).Cost = Patches(i).Scout(1).Cost;
    Patch(i).counter = Patches(i).Scout(1).counter;

end
%min_dist_idx = validx;

[~, recOrder]=sort(Bee_Recruitment,'descend');
Bee_Recruitment = Bee_Recruitment(recOrder);

size = linspace(1,1,min_dist_idx);

%% Bees Algorithm Local and Global Search
for it = 1:MaxIt
    if counter >= MaxEval
        break;
    end
    % All Sites (Exploitation and Exploration)
    for i = 1:min_dist_idx
        bestWorker.Cost = inf;
        assigntment = D_Triangular (0,size(i),1,1,Bee_Recruitment(i));
        for j = 1:Bee_Recruitment(i)
            [Pos1, Pos2, Pos3, Pos4] = Forage_Best(Patch(i).Position,assigntment(j)* Dims);
            Cost1 = ObjFunction(Pos1);
            Cost2 = ObjFunction(Pos2);
            Cost3 = ObjFunction(Pos3);
            Cost4 = ObjFunction(Pos4);
            [Worker.Position, Worker.Cost] = Find_Min(Pos1, Pos2, Pos3, Pos4, Cost1, Cost2, Cost3, Cost4);
            counter = counter+4;
            Worker.counter = counter;
            if Worker.Cost < bestWorker.Cost
                bestWorker = Worker;
            end
        end
        if bestWorker.Cost<Patch(i).Cost
            Patch(i) = bestWorker;
        end
    end

    % SORTING
    [~, RankOrder] = sort([Patch.Cost]);
    Patch = Patch(RankOrder);
    % Update Best Solution Ever Found
    OptSol = Patch(1);
    % taking of result
    OptCost(it) = OptSol.Cost;
    sonuc = OptSol.Cost;
    Counter(it) = counter;
    %Time(it) = toc;
    % Display Iteration Information
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(OptCost(it)) '; Fittness Evaluations = ' num2str(Counter(it))]);
    if(abs(Instance.optima-OptSol.Cost) == 0) 
        break;
    end
%     figure(1);
%     PlotSolution(OptSol.Position,Instance);
end
end