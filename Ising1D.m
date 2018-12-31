% This is a programme for one dimensional Ising model based on Metropolis algorithm which 
% is described by "Newman, M., & Barkema, G. (1999). Monte carlo methods in statistical 
% physics chapter 1-4. Oxford University Press: New York, USA."
% This algorithm is based on single-spin-flip dynamics.

t = 5000 ; % Time steps
N = 100 ; % The number of lattices 
T = 2.3 ; % The temperature of Ising model
k = 1 ; % Boltzmann¡¯s constant
J = 1 ; % Local interaction
B = 0 ; % Global effect
beta = 1/(k*T) ; 
model = struct ;
for i = 1:t
    model(i,1).state = zeros(1,N) ; % The states of Ising model, and the initial state (data(1,1).state) is set to condition when T=0.
    model(i,1).evolution = zeros(1,N+2) ; % The evolutions of Ising model. N+2 is for periodic boundary condition.
end
% Initial condition
model(1,1).state = model(1,1).state + 1 ;
model(1,1).evolution = model(1,1).evolution + 1 ;

for i = 1:t-1
    model(i+1,1).state = model(i,1).state ;
    model(i+1,1).evolution = model(i,1).evolution ;
    Index1 = fix(rand(1)*N + 1) ; % select a random location on the lattices
    West = model(i,1).evolution(Index1) ;
    East = model(i,1).evolution(Index1+2) ;
    Neighbour = [West,East] ;
    EnergyDiff = 2*J*model(i,1).state(Index1)*sum(Neighbour) + 2*B*model(i,1).state(Index1) ; 
    if EnergyDiff <= 0
        model(i+1,1).evolution(Index1+1) = -model(i,1).evolution(Index1+1) ;
    elseif EnergyDiff > 0
        Probability = exp(-beta*EnergyDiff) ;
        Seed = rand(1) ;
        if Seed < Probability
            model(i+1,1).evolution(Index1+1) = -model(i,1).evolution(Index1+1) ;
        elseif Seed >= Probability
            model(i+1,1).evolution(Index1+1) = model(i,1).evolution(Index1+1) ;
        end
    end
    model(i+1,1).state = model(i+1,1).evolution(2:N+1) ;
    % Periodic boundary condition
    model(i+1,1).evolution(1) = model(i+1,1).state(N) ;
    model(i+1,1).evolution(N+2) = model(i+1,1).state(1) ;
end

% Watch the results
A = zeros(t,N) ;
for i = 1:t
    A(i,:) = model(i,1).state ;
end
figure(2)
imagesc(A) ;