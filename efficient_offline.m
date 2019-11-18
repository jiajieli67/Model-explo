%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      +----------------------------+
%      |  Projection-based Model    |
%      |  online-efficient          |
%      +----------------------------+
% 
% Author: LI Jiajie 
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
