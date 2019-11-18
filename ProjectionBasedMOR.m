%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      +----------------------------+
%      |  Projection-based Model    |
%      |  order reduction methods   |
%      +----------------------------+
% 
% Author: Olivier Zahm (olivier.zahm@inria.fr)
% Date: November 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
set(groot,'defaulttextinterpreter','latex');  

% Initialize the model
model = parametrizedPDE();

%% Problem setup
model.plotPDEsetup()

%% Plot the solution
% Random draw of the parameter
X = model.randX();

% Compute the solution
u = model.u(X);

% Plot the solution
clf %->to clear the figure...
model.plotsol(u);

% Compute the quantity of interest:
Y = model.q' * u;


%% Quantity of interest as a random variable
% Compute K random draws of Y

K = 500;

Y = zeros(K,1);
tic
for k=1:K
    X = model.randX();
    Y(k) = model.q'*model.u(X);
end
time = toc;

disp('------------------------------')
disp('Time for computing x -> u(x) (full model) :')
disp(['    ' num2str(mean(time),3) ' sec per evaluation' ])


% Plot the histogram of Y
clf
hist(Y,30)
xlabel('Y')
ylabel('occurences of $Y$')


%% A (naive) construction of a reduced space
%%%%%%%%%%%%%%%%%
% OFFLINE PHASE %
%%%%%%%%%%%%%%%%%

% Span of r random snapshots of the solution:
r = 36;

tic; % start timer

% Compute r snapshots of x -> u(x)
Vr = zeros(model.n,r);
for k=1:r
    Vr(:,k) = model.u();
end
% Orthogonalize the columns of Vr:
Vr = orth(Vr);

timeOffline = toc; % stop timer

disp('------------------------------')
disp('Offine time (sec)')
disp(timeOffline)

%% Galerkin projection on the reduced space Vr
%%%%%%%%%%%%%%%%
% ONLINE PHASE %
%%%%%%%%%%%%%%%%

% Draw a new parameter value:
X = model.randX();

tic; % start timer

% Reduced operator and right-hand side
tildeA = Vr'*model.A(X)*Vr;
tildeb = Vr'*model.b;

% Solve the reduced linear system
lambda = tildeA\tildeb;

% Reconstruct the solution and compute the QoI
lambda = Vr*lambda;
tildeY = model.q'*lambda;

timeOnline = toc; % stop timer

disp('------------------------------')
disp('Online time (sec)')
disp(timeOnline)

%% Compare with the full model

tic; % start timer

% Solve the full system and compute YX
u = model.A(X)\model.b;
Y = model.q'*u;

timeFullModel = toc; % stop timer

disp('------------------------------')
disp('Time for evaluating the full model (sec)')
disp(timeFullModel)

disp(' Speed-up:')
disp(timeFullModel/timeOnline)

disp('Error')
error = abs(Y - tildeY)/abs(Y);
disp(error)

% Plot the solutions
clf
subplot(1,2,1)
model.plotsol(u)
title('Full model')
subplot(1,2,2)
model.plotsol(lambda)
title('Reduced model')

%% Error over the parameter domain
% Compute the (exact) error on K random parameters

K = 500;

error = zeros(K,1);
for k=1:K
    % Take a new parameter value:
    X = model.randX();
    
    % Reduced solution
    tildeA = Vr'*model.A(X)*Vr;
    tildeb  = Vr'*model.b;
    lambda = tildeA\tildeb;
    utilde = Vr*lambda;
    
    % Full order solution
    u = model.u(X);
    
    % Evaluate the true error
    error(k) = norm(u-utilde);
end

% Plot the histograms of the error and of the residual
clf
hist( log10(error) ,15)
ylabel('occurence')
xlabel('error measure ($\log_{10}$-scale)')

disp('------------------------------')
disp('L-infty error:')
disp( max(error) )
disp('L-2 error:')
disp( sqrt( sum(error.^2)/K ) )

