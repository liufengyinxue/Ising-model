% This is a programme for one dimensional Ising model based on Metropolis algorithm which 
% is described by "Newman, M., & Barkema, G. (1999). Monte carlo methods in statistical 
% physics chapter 1-4. Oxford University Press: New York, USA."
% This algorithm is based on single-spin-flip dynamics.

% This programme samples the results of 2D Ising model with the interval 'dt' 

t = 1e+7 ; % Time steps
M = 100 ; % The row of lattice
N = 100 ; % The column of lattice 
T = 2.4 ; % The temperature of Ising model
k = 1 ; % Boltzmann¡¯s constant
J = 1 ; % Local interaction
B = 0 ; % Global effect
beta = 1/(k*T) ;
dt = 1e+3 ; % Sample interval
Burntime = 1e+5 ; % Non-sampling time
SamNum = fix((t - Burntime)/dt) ; % Numnber of samples
% Save samples of result
model = struct ;
for j = 1:SamNum
    model(j,1).data = zeros(M,N) ; % Samples of result
end
% Initial condition
modelstate = zeros(M,N) + 1 ;
modelevolution = zeros(M+2,N+2) + 1 ;

for i = 1:t
    % Select a random location on the lattices
    Index1 = fix(rand(1)*M + 1) ;
    Index2 = fix(rand(1)*N + 1) ;
    % Find neighbours
    North = modelevolution(Index1,Index2+1) ;
    South = modelevolution(Index1+2,Index2+1) ;
    East = modelevolution(Index1+1,Index2+2) ;
    West = modelevolution(Index1+1,Index2) ;
    Neighbour = [North,South,East,West] ;
    % Calculate energy difference
    EnergyDiff = 2*J*modelstate(Index1,Index2)*sum(Neighbour) + 2*B*modelstate(Index1,Index2) ; 
    % Flip or not flip
    if EnergyDiff <= 0
        modelevolution(Index1+1,Index2+1) = -modelevolution(Index1+1,Index2+1) ;
    elseif EnergyDiff > 0
        Probability = exp(-beta*EnergyDiff) ;
        Seed = rand(1) ;
        if Seed < Probability
            modelevolution(Index1+1,Index2+1) = -modelevolution(Index1+1,Index2+1) ;
        elseif Seed >= Probability
            modelevolution(Index1+1,Index2+1) = modelevolution(Index1+1,Index2+1) ;
        end
    end
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

% Watch the results
for j = 1:SamNum
    A = model(j,1).data ;
    imagesc(A) ;
    M = getframe;
end
movie(M,1,10) %Play once at 10 frames per second