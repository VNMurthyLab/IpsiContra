%% makeCrossCortex.m
% Program to create a cross cortical linking connectivity matrix.
%
% INPUT
% alpha - moves between random and structured, 0 purely random, 1 purely
% structured, 0.5 half and half.
% yA (ipsi) and yB (contra) are the two representations on either side of the cortex
% th is the threshold for the neurons
% G is the connectivity matrix, in case it has already been calculated
%
% OUTPUT
% z is the represention of the contra on the ipsi side
% GRand is the random connectivity matrix
% GStruct is the structured connectivity matrix

%%
function [z, GRand, GStruct] = makeCrossCortex(alpha, yA, yB, th, G)
global Ny Gc SContra

NumberOdours = size(yA, 2);

if ~exist('th', 'var') || isempty(th)
    error('Should give threshold!')
end

makeG = false;
if ~exist('G','var') || isempty(G)
    makeG = true;
end

if makeG
    if alpha > 0
        % Here we will start doing work to make the correlated G matrix
        % First we need the two C vectors, the average over the odours of the reps
        CA = 1/NumberOdours*sum(yA, 2);
        CB = 1/NumberOdours*sum(yB, 2);
        
        % Now we sum up the correlated parts
        Summing = zeros(Ny, Ny);
        for i = 1:NumberOdours
            Summing = Summing + (yA(:, i)-CA)*(yB(:, i)-CB)';
        end
        GStruct = 1/(NumberOdours)*Summing;%*1/Ny
        clear Summing
    end
    
    if alpha ~= 1
        % Now we start making the random G matrix
        % Lets just say it is a gaussian random matrix with a given sparseness
        GLoc = (rand(Ny, Ny)< Gc);
        GRand = random('normal', 0, 1, [Ny, Ny]);
        GRand(~GLoc) = 0;
    end
     
    % Before finally constructing the matrix
    if alpha == 0
        G = GRand;
    elseif alpha == 1
        G = GStruct;
    else
        % Here we combine them but additionally normalise such that the
        % magnitude of the outputs of the two matrices are the same
        StructNorm = norm(GStruct*yB, 2);
        RandNorm = norm(GRand*yB, 2);
        psi = RandNorm/StructNorm;
        GStruct = GStruct*psi;
        G = alpha*GStruct + (1-alpha)*GRand;
    end
end

% Now propagate the contra rep through the matrix
z = G*yB;

% Now we find out how much we have to rescale the representations in order
% that the inputted threshold achieves the desired sparseness (SContra)
sz = size(z);
No = sz(2);
temp = sort(reshape(abs(z),sz(1)*sz(2)*size(z,3),1),'descend');
loc = round(Ny*SContra*No); % Location of the last active neuron we want to choose active
phi = th/temp(loc);
% An error I have never had the misfortune of having to implement a
% correction for.
if temp(loc) < 0 || temp(loc) == 0
    disp('ERROR! (OF SORTS): the threshold activity is negative, requires perhaps more thought to make this work')
end

% Now rescale all the other things as desired
if alpha ~= 1 && makeG
    GRand = phi*GRand;
end
if alpha > 0 && makeG
    GStruct = phi*GStruct;
end
z = phi*z;

% Before finally thresholding this representation
z1 = z-th;
z1(z1<0)=0;
z2 = z+th;
z2(z2>0)=0;
z = z1 + z2;
end