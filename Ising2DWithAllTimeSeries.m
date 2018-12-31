% This is a programme for one dimensional Ising model based on Metropolis algorithm which 
% is described by "Newman, M., & Barkema, G. (1999). Monte carlo methods in statistical 
% physics chapter 1-4. Oxford University Press: New York, USA."
% This algorithm is based on single-spin-flip dynamics.

% Note:not recommend use this programme because the computers in 2018 don't
% have the ability to save so much data. Please go to "Ising2D.m" to
% execute 2D Ising model.

t = 1e+4 ; % Time steps
M = 100 ; % The row of lattice
N = 100 ; % The column of lattice 
T = 2.4 ; % The temperature of Ising model
k = 1 ; % Boltzmann’s constant
J = 1 ; % Local interaction
B = 0 ; % Global effect
beta = 1/(k*T) ; 
model = struct ;
for i = 1:t
    model(i,1).state = zeros(M,N) ; % The states of Ising model, and the initial state (data(1,1).state) is set to condition when T=0.
    model(i,1).evolution = zeros(M+2,N+2) ; % The evolutions of Ising model. N+2 is for periodic boundary condition.
end
% Initial condition
model(1,1).state = model(1,1).state + 1 ;
model(1,1).evolution = model(1,1).evolution + 1 ;

for i = 1:t-1
    model(i+1,1).state = model(i,1).state ;
    model(i+1,1).evolution = model(i,1).evolution ;
    % Select a random location on the lattices
    Index1 = fix(rand(1)*M + 1) ;
    Index2 = fix(rand(1)*N + 1) ;
    % Find neighbours
    North = model(i,1).evolution(Index1,Index2+1) ;
    South = model(i,1).evolution(Index1+2,Index2+1) ;
    East = model(i,1).evolution(Index1+1,Index2+2) ;
    West = model(i,1).evolution(Index1+1,Index2) ;
    Neighbour = [North,South,East,West] ;
    % Calculate energy difference
    EnergyDiff = 2*J*model(i,1).state(Index1,Index2)*sum(Neighbour) + 2*B*model(i,1).state(Index1,Index2) ; 
    % Flip or not flip
    if EnergyDiff <= 0
        model(i+1,1).evolution(Index1+1,Index2+1) = -model(i,1).evolution(Index1+1,Index2+1) ;
    elseif EnergyDiff > 0
        Probability = exp(-beta*EnergyDiff) ;
        Seed = rand(1) ;
        if Seed < Probability
            model(i+1,1).evolution(Index1+1,Index2+1) = -model(i,1).evolution(Index1+1,Index2+1) ;
        elseif Seed >= Probability
            model(i+1,1).evolution(Index1+1,Index2+1) = model(i,1).evolution(Index1+1,Index2+1) ;
        end
    end
    % Save state from evolution
    model(i+1,1).state = model(i+1,1).evolution(2:M+1,2:N+1) ;
    % Periodic boundary condition
    model(i+1,1).evolution(1,2:N+1) = model(i+1,1).state(M,:) ;  
    model(i+1,1).evolution(M+2,2:N+1) = model(i+1,1).state(1,:) ;  
    model(i+1,1).evolution(2:N+1,1) = model(i+1,1).state(:,N) ;  
    model(i+1,1).evolution(2:N+1,N+2) = model(i+1,1).state(:,1) ;
    model(i+1,1).evolution(1,1) = model(i+1,1).state(M,N) ; 
    model(i+1,1).evolution(1,N+2) = model(i+1,1).state(M,1) ; 
    model(i+1,1).evolution(M+2,1) = model(i+1,1).state(1,N) ; 
    model(i+1,1).evolution(M+2,N+2) = model(i+1,1).state(1,1) ;
end

% % Watch the results
% for j = 1:t
%     A = model(j,1).state ;
%     imagesc(A) ;
%     M = getframe;
% end
% movie(M,1,10) %每秒10帧的速度播放1次