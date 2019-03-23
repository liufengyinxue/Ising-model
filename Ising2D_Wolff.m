% This is a programme for two dimensional Ising model based on Wolff algorithm which 
% is described by "Newman, M., & Barkema, G. (1999). Monte carlo methods in statistical 
% physics. Oxford University Press: New York, USA."
% This algorithm is based on cluster-flip algorithm.
% This programme samples the results of 2D Ising model with the interval 'dt' 

t = 1e+7 ; % Time steps
M = 100 ; % The row of lattice
N = 100 ; % The column of lattice 
T = 2.3 ; % The temperature of Ising model
k = 1 ; % Boltzmann¡¯s constant
J = 1 ; % Local interaction
B = 0 ; % Global effect
beta = 1/(k*T) ;
Padd = 1 - exp(-2*beta*J) ;
dt = 1e+3 ; % Sample interval
Burntime = 1e+5 ; % Non-sampling time
SamNum = fix((t - Burntime)/dt) ; % Number of samples
model = struct ;% Save samples of result
for j = 1:SamNum
    model(j,1).data = zeros(M,N) ; % Samples of result
end
% Initial condition
modelstate = zeros(M,N) + 1 ;
modelevolution = zeros(M+2,N+2) + 1 ;

for i = 1:t
    % Select a random seed in the lattice
    Index1 = fix(rand(1)*M + 1) ;
    Index2 = fix(rand(1)*N + 1) ;
    Seed = [Index1,Index2] ;
    Stack = [] ; % Save index for model.state
    Stack = cat(1,Stack,Seed) ;
    Cluster = [] ; % Save index for model.state
    Cluster = cat(1,Cluster,Seed) ;
    while ~isempty(Stack)
        [m1,n1] = size(Stack) ;
        for j = 1:m1
            if Stack(j,1) == 0 && Stack(j,2) >= 1 && Stack(j,2) <= N
                Centre = modelevolution(Stack(j,1)+1,Stack(j,2)+1) ;
                South = modelevolution(Stack(j,1)+2,Stack(j,2)+1) ;
                SN = (Stack(j,1) + 1 - Cluster(:,1)).^2 + (Stack(j,2) - Cluster(:,2)).^2 ;
                S1 = length(find(SN == 0)) ;
                East = modelevolution(Stack(j,1)+1,Stack(j,2)+2) ;
                EN = (Stack(j,1) - Cluster(:,1)).^2 + (Stack(j,2) + 1 - Cluster(:,2)).^2 ;
                E1 = length(find(EN == 0)) ;
                West = modelevolution(Stack(j,1)+1,Stack(j,2)) ;
                WN = (Stack(j,1) - Cluster(:,1)).^2 + (Stack(j,2) - 1 - Cluster(:,2)).^2 ;
                W1 = length(find(WN == 0)) ;
                Neighbour = [0,South,East,West] ;
                CeNe = Centre*Neighbour ;
                S = rand(4,1) ;
                if CeNe(2) == 1 && S(2) <= Padd && S1 == 0
                    tempS = [Stack(j,1)+1,Stack(j,2)] ;
                    Stack = cat(1,Stack,tempS) ;
                    Cluster = cat(1,Cluster,tempS) ;
                end
                if CeNe(3) == 1 && S(3) <= Padd && E1 == 0
                    tempE = [Stack(j,1),Stack(j,2)+1] ;
                    Stack = cat(1,Stack,tempE) ;
                    Cluster = cat(1,Cluster,tempE) ;
                end
                if CeNe(4) == 1 && S(4) <= Padd && W1 == 0
                    tempW = [Stack(j,1),Stack(j,2)-1] ;
                    Stack = cat(1,Stack,tempW) ;
                    Cluster = cat(1,Cluster,tempW) ;
                end
            elseif Stack(j,1) == M+1 && Stack(j,2) >= 1 && Stack(j,2) <= N
                Centre = modelevolution(Stack(j,1)+1,Stack(j,2)+1) ;
                North = modelevolution(Stack(j,1),Stack(j,2)+1) ;
                NN = (Stack(j,1) - 1 - Cluster(:,1)).^2 + (Stack(j,2) - Cluster(:,2)).^2 ;
                N1 = length(find(NN == 0)) ;
                East = modelevolution(Stack(j,1)+1,Stack(j,2)+2) ;
                EN = (Stack(j,1) - Cluster(:,1)).^2 + (Stack(j,2) + 1 - Cluster(:,2)).^2 ;
                E1 = length(find(EN == 0)) ;
                West = modelevolution(Stack(j,1)+1,Stack(j,2)) ;
                WN = (Stack(j,1) - Cluster(:,1)).^2 + (Stack(j,2) - 1 - Cluster(:,2)).^2 ;
                W1 = length(find(WN == 0)) ;
                Neighbour = [North,0,East,West] ;
                CeNe = Centre*Neighbour ;
                S = rand(4,1) ;
                if CeNe(1) == 1 && S(1) <= Padd && N1 == 0 
                    tempN = [Stack(j,1) - 1,Stack(j,2)] ;
                    Stack = cat(1,Stack,tempN) ;
                    Cluster = cat(1,Cluster,tempN) ;
                end
                if CeNe(3) == 1 && S(3) <= Padd && E1 == 0
                    tempE = [Stack(j,1),Stack(j,2)+1] ;
                    Stack = cat(1,Stack,tempE) ;
                    Cluster = cat(1,Cluster,tempE) ;
                end
                if CeNe(4) == 1 && S(4) <= Padd && W1 == 0
                    tempW = [Stack(j,1),Stack(j,2)-1] ;
                    Stack = cat(1,Stack,tempW) ;
                    Cluster = cat(1,Cluster,tempW) ;
                end        
            elseif Stack(j,2) == 0 && Stack(j,1) >= 1 && Stack(j,1) <= M
                Centre = modelevolution(Stack(j,1)+1,Stack(j,2)+1) ;
                North = modelevolution(Stack(j,1),Stack(j,2)+1) ;
                NN = (Stack(j,1) - 1 - Cluster(:,1)).^2 + (Stack(j,2) - Cluster(:,2)).^2 ;
                N1 = length(find(NN == 0)) ;
                South = modelevolution(Stack(j,1)+2,Stack(j,2)+1) ;
                SN = (Stack(j,1) + 1 - Cluster(:,1)).^2 + (Stack(j,2) - Cluster(:,2)).^2 ;
                S1 = length(find(SN == 0)) ;
                East = modelevolution(Stack(j,1)+1,Stack(j,2)+2) ;
                EN = (Stack(j,1) - Cluster(:,1)).^2 + (Stack(j,2) + 1 - Cluster(:,2)).^2 ;
                E1 = length(find(EN == 0)) ;
                Neighbour = [North,South,East,0] ;
                CeNe = Centre*Neighbour ;
                S = rand(4,1) ;
                if CeNe(1) == 1 && S(1) <= Padd && N1 == 0 
                    tempN = [Stack(j,1) - 1,Stack(j,2)] ;
                    Stack = cat(1,Stack,tempN) ;
                    Cluster = cat(1,Cluster,tempN) ;
                end
                if CeNe(2) == 1 && S(2) <= Padd && S1 == 0
                    tempS = [Stack(j,1)+1,Stack(j,2)] ;
                    Stack = cat(1,Stack,tempS) ;
                    Cluster = cat(1,Cluster,tempS) ;
                end
                if CeNe(3) == 1 && S(3) <= Padd && E1 == 0
                    tempE = [Stack(j,1),Stack(j,2)+1] ;
                    Stack = cat(1,Stack,tempE) ;
                    Cluster = cat(1,Cluster,tempE) ;
                end       
            elseif Stack(j,2) == N+1 && Stack(j,1) >= 1 && Stack(j,1) <= M
                Centre = modelevolution(Stack(j,1)+1,Stack(j,2)+1) ;
                North = modelevolution(Stack(j,1),Stack(j,2)+1) ;
                NN = (Stack(j,1) - 1 - Cluster(:,1)).^2 + (Stack(j,2) - Cluster(:,2)).^2 ;
                N1 = length(find(NN == 0)) ;
                South = modelevolution(Stack(j,1)+2,Stack(j,2)+1) ;
                SN = (Stack(j,1) + 1 - Cluster(:,1)).^2 + (Stack(j,2) - Cluster(:,2)).^2 ;
                S1 = length(find(SN == 0)) ;               
                West = modelevolution(Stack(j,1)+1,Stack(j,2)) ;
                WN = (Stack(j,1) - Cluster(:,1)).^2 + (Stack(j,2) - 1 - Cluster(:,2)).^2 ;
                W1 = length(find(WN == 0)) ;
                Neighbour = [North,South,0,West] ;
                CeNe = Centre*Neighbour ;
                S = rand(4,1) ;
                if CeNe(1) == 1 && S(1) <= Padd && N1 == 0 
                    tempN = [Stack(j,1) - 1,Stack(j,2)] ;
                    Stack = cat(1,Stack,tempN) ;
                    Cluster = cat(1,Cluster,tempN) ;
                end
                if CeNe(2) == 1 && S(2) <= Padd && S1 == 0
                    tempS = [Stack(j,1)+1,Stack(j,2)] ;
                    Stack = cat(1,Stack,tempS) ;
                    Cluster = cat(1,Cluster,tempS) ;
                end
                if CeNe(4) == 1 && S(4) <= Padd && W1 == 0
                    tempW = [Stack(j,1),Stack(j,2)-1] ;
                    Stack = cat(1,Stack,tempW) ;
                    Cluster = cat(1,Cluster,tempW) ;
                end            
            elseif Stack(j,1) == 0 && Stack(j,2) == 0 
                 Centre = modelevolution(Stack(j,1)+1,Stack(j,2)+1) ;
                 South = modelevolution(Stack(j,1)+2,Stack(j,2)+1) ;
                 SN = (Stack(j,1) + 1 - Cluster(:,1)).^2 + (Stack(j,2) - Cluster(:,2)).^2 ;
                 S1 = length(find(SN == 0)) ;
                 East = modelevolution(Stack(j,1)+1,Stack(j,2)+2) ;
                 EN = (Stack(j,1) - Cluster(:,1)).^2 + (Stack(j,2) + 1 - Cluster(:,2)).^2 ;
                 E1 = length(find(EN == 0)) ;
                 Neighbour = [0,South,East,0] ;
                 CeNe = Centre*Neighbour ;
                 S = rand(4,1) ;
                 if CeNe(2) == 1 && S(2) <= Padd && S1 == 0
                     tempS = [Stack(j,1)+1,Stack(j,2)] ;
                     Stack = cat(1,Stack,tempS) ;
                     Cluster = cat(1,Cluster,tempS) ;
                 end
                 if CeNe(3) == 1 && S(3) <= Padd && E1 == 0
                     tempE = [Stack(j,1),Stack(j,2)+1] ;
                     Stack = cat(1,Stack,tempE) ;
                     Cluster = cat(1,Cluster,tempE) ;
                 end        
            elseif Stack(j,1) == 0 && Stack(j,2) == N+1
                 Centre = modelevolution(Stack(j,1)+1,Stack(j,2)+1) ;
                 South = modelevolution(Stack(j,1)+2,Stack(j,2)+1) ;
                 SN = (Stack(j,1) + 1 - Cluster(:,1)).^2 + (Stack(j,2) - Cluster(:,2)).^2 ;
                 S1 = length(find(SN == 0)) ;
                 West = modelevolution(Stack(j,1)+1,Stack(j,2)) ;
                 WN = (Stack(j,1) - Cluster(:,1)).^2 + (Stack(j,2) - 1 - Cluster(:,2)).^2 ;
                 W1 = length(find(WN == 0)) ;
                 Neighbour = [0,South,0,West] ;
                 CeNe = Centre*Neighbour ;
                 S = rand(4,1) ;
                 if CeNe(2) == 1 && S(2) <= Padd && S1 == 0
                     tempS = [Stack(j,1)+1,Stack(j,2)] ;
                     Stack = cat(1,Stack,tempS) ;
                     Cluster = cat(1,Cluster,tempS) ;
                 end
                 if CeNe(4) == 1 && S(4) <= Padd && W1 == 0
                     tempW = [Stack(j,1),Stack(j,2)-1] ;
                     Stack = cat(1,Stack,tempW) ;
                     Cluster = cat(1,Cluster,tempW) ;
                 end 
            elseif Stack(j,1) == M+1 && Stack(j,2) == 0
                Centre = modelevolution(Stack(j,1)+1,Stack(j,2)+1) ;
                North = modelevolution(Stack(j,1),Stack(j,2)+1) ;
                NN = (Stack(j,1) - 1 - Cluster(:,1)).^2 + (Stack(j,2) - Cluster(:,2)).^2 ;
                N1 = length(find(NN == 0)) ;
                East = modelevolution(Stack(j,1)+1,Stack(j,2)+2) ;
                EN = (Stack(j,1) - Cluster(:,1)).^2 + (Stack(j,2) + 1 - Cluster(:,2)).^2 ;
                E1 = length(find(EN == 0)) ;
                Neighbour = [North,0,East,0] ;
                CeNe = Centre*Neighbour ;
                S = rand(4,1) ;
                if CeNe(1) == 1 && S(1) <= Padd && N1 == 0 
                    tempN = [Stack(j,1) - 1,Stack(j,2)] ;
                    Stack = cat(1,Stack,tempN) ;
                    Cluster = cat(1,Cluster,tempN) ;
                end
                if CeNe(3) == 1 && S(3) <= Padd && E1 == 0
                    tempE = [Stack(j,1),Stack(j,2)+1] ;
                    Stack = cat(1,Stack,tempE) ;
                    Cluster = cat(1,Cluster,tempE) ;
                end
            elseif Stack(j,1) == M+1 && Stack(j,2) == N+1
                Centre = modelevolution(Stack(j,1)+1,Stack(j,2)+1) ;
                North = modelevolution(Stack(j,1),Stack(j,2)+1) ;
                NN = (Stack(j,1) - 1 - Cluster(:,1)).^2 + (Stack(j,2) - Cluster(:,2)).^2 ;
                N1 = length(find(NN == 0)) ;
                West = modelevolution(Stack(j,1)+1,Stack(j,2)) ;
                WN = (Stack(j,1) - Cluster(:,1)).^2 + (Stack(j,2) - 1 - Cluster(:,2)).^2 ;
                W1 = length(find(WN == 0)) ;
                Neighbour = [North,0,0,West] ;
                CeNe = Centre*Neighbour ;
                S = rand(4,1) ;
                if CeNe(1) == 1 && S(1) <= Padd && N1 == 0 
                    tempN = [Stack(j,1) - 1,Stack(j,2)] ;
                    Stack = cat(1,Stack,tempN) ;
                    Cluster = cat(1,Cluster,tempN) ;
                end 
                if CeNe(4) == 1 && S(4) <= Padd && W1 == 0
                    tempW = [Stack(j,1),Stack(j,2)-1] ;
                    Stack = cat(1,Stack,tempW) ;
                    Cluster = cat(1,Cluster,tempW) ;
                end            
            elseif Stack(j,1) >= 1 && Stack(j,1) <= M && Stack(j,2) >= 1 && Stack(j,2) <= N
                Centre = modelevolution(Stack(j,1)+1,Stack(j,2)+1) ;
                North = modelevolution(Stack(j,1),Stack(j,2)+1) ;
                NN = (Stack(j,1) - 1 - Cluster(:,1)).^2 + (Stack(j,2) - Cluster(:,2)).^2 ;
                N1 = length(find(NN == 0)) ;
                South = modelevolution(Stack(j,1)+2,Stack(j,2)+1) ;
                SN = (Stack(j,1) + 1 - Cluster(:,1)).^2 + (Stack(j,2) - Cluster(:,2)).^2 ;
                S1 = length(find(SN == 0)) ;
                East = modelevolution(Stack(j,1)+1,Stack(j,2)+2) ;
                EN = (Stack(j,1) - Cluster(:,1)).^2 + (Stack(j,2) + 1 - Cluster(:,2)).^2 ;
                E1 = length(find(EN == 0)) ;
                West = modelevolution(Stack(j,1)+1,Stack(j,2)) ;
                WN = (Stack(j,1) - Cluster(:,1)).^2 + (Stack(j,2) - 1 - Cluster(:,2)).^2 ;
                W1 = length(find(WN == 0)) ;
                Neighbour = [North,South,East,West] ;
                CeNe = Centre*Neighbour ;
                S = rand(4,1) ;
                if CeNe(1) == 1 && S(1) <= Padd && N1 == 0 
                    tempN = [Stack(j,1) - 1,Stack(j,2)] ;
                    Stack = cat(1,Stack,tempN) ;
                    Cluster = cat(1,Cluster,tempN) ;
                end
                if CeNe(2) == 1 && S(2) <= Padd && S1 == 0
                    tempS = [Stack(j,1)+1,Stack(j,2)] ;
                    Stack = cat(1,Stack,tempS) ;
                    Cluster = cat(1,Cluster,tempS) ;
                end
                if CeNe(3) == 1 && S(3) <= Padd && E1 == 0
                    tempE = [Stack(j,1),Stack(j,2)+1] ;
                    Stack = cat(1,Stack,tempE) ;
                    Cluster = cat(1,Cluster,tempE) ;
                end
                if CeNe(4) == 1 && S(4) <= Padd && W1 == 0
                    tempW = [Stack(j,1),Stack(j,2)-1] ;
                    Stack = cat(1,Stack,tempW) ;
                    Cluster = cat(1,Cluster,tempW) ;
                end  
            end          
        end
        Stack(1:m1,:) = [] ;
    end
    modelevolution(Cluster(:,1)+1,Cluster(:,2)+1) = -modelevolution(Cluster(:,1)+1,Cluster(:,2)+1) ; 
   % Save state from evolution
    modelstate = modelevolution(2:M+1,2:N+1) ;
   % Sample
    if (i - Burntime) > 0 && rem((i - Burntime)/dt,1) == 0
        model((i - Burntime)/dt,1).data = modelstate ;
    end
    % Periodic boundary condition
    modelevolution(1,2:N+1) = modelstate(M,:) ;  
    modelevolution(M+2,2:N+1) = modelstate(1,:) ;  
    modelevolution(2:N+1,1) = modelstate(:,N) ;  
    modelevolution(2:N+1,N+2) = modelstate(:,1) ;
    modelevolution(1,1) = modelstate(M,N) ; 
    modelevolution(1,N+2) = modelstate(M,1) ; 
    modelevolution(M+2,1) = modelstate(1,N) ; 
    modelevolution(M+2,N+2) = modelstate(1,1) ; 
end

% % Watch the results
% for j = 1:SamNum
%     A = model(j,1).data ;
%     imagesc(A) ;
%     M = getframe;
% end
% movie(M,1,10) %Play once at 10 frames per second


