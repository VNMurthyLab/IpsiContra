%% PiriformModelling.m
%
% In this file we run simulations of the piriform and use them to compare
% to the data. With this we attempt to validate the idea that piriform
% cortical representations could be aligned by a structured cross-cortical
% matrix.
%
% The structure of the file is as follows:
% 1) Load things
% 2) Create the representations for ipsi and contra with varying
%       cross-cortical connectivity.
% 3) Extract four measures of alignment between the representations:
%   i) Percentage of bilaterally correlated neurons
%   ii) Correlation of ipsi and contra responses
%   iii) Ipsi decoding contra odour ID accuracy
%   iv) Side decoding accuracy
% 4) Then compare to data which allows us to pick the cross-cortical
%       connectivity that fits best
% 5) Before plotting a lot.

%% Load and setup variables
% Setup variables used throughout the simulations:
%   SIpsi, SContra - Sparseness of Ipsi and Contra representations
%   Ny - number of neurons to simulate in each piriform cortex
%   Gc - sparseness of cross-cortical random connectivity
global SIpsi SContra Ny  Gc
load('Data_and_Notes\NeuronalDataforModelComparison.mat')
Ny = 5000;
SIpsi = S(2);
SContra = S(1);
Gc = 0.1;
NumNeuExp = 385;
Odours = 15;

% The column convention between the loaded data and the simulations is,
% unfortunately different, so lets change that quickly.
AON = permute(AON, [2, 4, 1, 3]);

%% Now create piriform representations
% Get correlated multivariate normal, one for each side
h1 = mvnrnd(zeros(Odours,1), CosineSim, Ny);
h2 = mvnrnd(zeros(Odours,1), CosineSim, Ny);

% Threshold the shit out of that bad boy
phi = sqrt(2)*erfcinv(SIpsi);
h1(abs(h1) < phi) = 0;
h2(abs(h2) < phi) = 0;
h1(abs(h1)>0) = h1(abs(h1)>0) - sign(h1(abs(h1)>0)).*phi;
h2(abs(h2)>0) = h2(abs(h2)>0) - sign(h2(abs(h2)>0)).*phi; 

% Now for a selection of alphas we create the contra representation
alphas = [0:0.05:0.35, 0.36:0.01:0.44, 0.45:0.05:1];
z = zeros(Ny,Odours,length(alphas));
counter = 0; 
% counts through options:
% 0 means you are using only GRand therefore no need for special effort
% 1 means you need to make the structured cross-cortical matrix and store it!
% 2 means the structured matrix has already been made, and you should use
%       use it to make the cross-cortical matrix (G) and pass as input.
for i = 1:length(alphas)
    alpha = alphas(i);
    disp(['Alpha: ', num2str(alpha)])
    if alpha > 0 && counter < 2
        counter = counter +  1;
    end
    if counter == 0
        z(:,:,i) = makeCrossCortex(alpha, h1, h2, phi);
    elseif counter == 1
        [z(:,:,i), GRand, GStruct] = makeCrossCortex(alpha, h1, h2, phi);
    else
        G = (1-alpha)*GRand + alpha*GStruct;
        z(:,:,i) = makeCrossCortex(alpha, h1, h2, phi, G);
    end
end

% Finally we correct the units to firing rate
h1 = h1/ScaleConversion;
z = z/ScaleConversion;

%% Here we measure the % bilaterally correlated neurons in the population
% We will subsample 385 neurons (number measured in experiment) and
% calculated the % of those that are bilaterally correlated. Num controls
% how many of these subsampling we do.
%
% THIS ONE TAKES A WHILE!
Num = 200;
BilCorr = zeros(Num, length(alphas));       % Will hold %s

for j = 1:length(alphas)
    for k = 1:Num
        % Choose a random set of neurons
        set = randperm(Ny, NumNeuExp);
        Contra = z(set,:,j);
        Ipsi = h1(set,:);
        counter = 0;
        for i = 1:385
            NeuronIpsi = Ipsi(i,:);
            NeuronContra = Contra(i,:);
            [R,P] = corrcoef(NeuronIpsi, NeuronContra);
            if P(1,2) < 0.05
                counter = counter + 1;
            end
        end
        BilCorr(k,j) = counter/385;
    end
end

BilCorrDev = std(BilCorr, 0, 1);
BilCorr = mean(BilCorr, 1);

% Elsewhere in the project it was calculated that 135 neurons in AON were
% bilaterally correlated.
BilCorrData = 135/NumNeuExp;

%% Second measurement, the correlation coefficient between ipsi and contra
% Again across Num iterations of a 385 neuron subsampling
CorrCoef = zeros(Num, length(alphas));

for mu = 1:length(alphas)
    for j = 1:Num
        set = randperm(Ny, NumNeuExp);
        resp = z(set,:,mu);
        resp2 = h1(set, :);
        CorrCoef(j,mu) = corr(resp2(:), resp(:));
    end
end

CorrCoefDev = std(CorrCoef, 0, 1);
CorrCoef = mean(CorrCoef, 1);

% Same for data
CorrCoefData = corr(reshape(AON(:,:,:,1), [floor(numel(AON)/2),1]),reshape(AON(:,:,:,2), [floor(numel(AON)/2),1]));

%% Third we decode contra odour ID using ipsi and vice versa for model and data.
% We use a different classification technique here vs elsewhere.
% This is because the model does not have noise/a trial structure.
% We therefore simply classify odour centroids from one side to their
% nearest neighbour from the other and ask for % correctly assigned.
NumNeus = [0:10:50, 100:50:NumNeuExp, NumNeuExp];   % Neurons to use
Indices = length(NumNeus);
Accuracy = zeros(length(alphas), length(NumNeus),2);
Deviation = Accuracy;

% This method is implemented in CentroidClassifier which we call here.
for i = 1:length(alphas)
    resp = squeeze(z(:,:,i));
    for k = 1:length(NumNeus)
        NumNeu = NumNeus(k);
        [Accuracy(i,k,1), Deviation(i,k,1)] = CentroidClassifier(h1, resp, NumNeu, Num);
        [Accuracy(i,k,2), Deviation(i,k,2)] = CentroidClassifier(resp, h1, NumNeu, Num);
    end
end

% Now we find the data centroids and perform the same task
IpsiCent = findCentroids(squeeze(AON(:,:,:,2)));
ContraCent = findCentroids(squeeze(AON(:,:,:,1)));
AccuracyData = zeros(length(NumNeus), 2);
DeviationData = AccuracyData;

for k = 1:length(NumNeus)
    NumNeu = NumNeus(k);
    [AccuracyData(k, 1), DeviationData(k, 1)] = CentroidClassifier(IpsiCent', ContraCent', NumNeu, Num);
    [AccuracyData(k, 2), DeviationData(k, 2)] = CentroidClassifier(ContraCent', IpsiCent', NumNeu, Num);
end

% For reasons we discuss a little in figure S7 sometimes angular nearest
% neighbour rather than euclidean is best. Here we implement that.
% (It is basically due to the mismatch in absolute size in contra vs ipsi
% representations that can sometimes result)
AccuracyAng = zeros(length(alphas),1);
DeviationAng = AccuracyAng;
[AccuracyData2, DeviationData2] = CentroidClassifier2(ContraCent', IpsiCent', 385, Num);
for i = 1:length(alphas)
    resp = squeeze(z(:,:,i));
    [AccuracyAng(i), DeviationAng(i)] = CentroidClassifier2(resp, h1, 385, Num);
end

%% Finally we try and decode the side of presentation
% Again we use our centroid classifying technique, but in this case it is
% just a binary comparison: is the query point closer to contra or ipsi centroid?
AccuracySide = zeros(length(alphas), length(NumNeus));
DeviationSide = AccuracySide;

for i = 1:length(alphas)
    resp = squeeze(z(:,:,i));
    for k = 1:length(NumNeus)
        NumNeu = NumNeus(k);
        [AccuracySide(i,k), DeviationSide(i,k)] = CentroidClassifySide(h1, resp, NumNeu, Num);
    end
end

% And calculate the same for the data
AccuracyDataSide = zeros(length(NumNeus),1);
DeviationDataSide = AccuracyDataSide;

for k = 1:length(NumNeus)
    NumNeu = NumNeus(k);
    [AccuracyDataSide(k), DeviationDataSide(k)] = CentroidClassifySide(IpsiCent', ContraCent', NumNeu, Num);
end    

%% Finally, for each alpha value find z-score and choose lowest
% The Z score is the ratio of distance between point and expected value
% squared over the standard deviation. We choose the alpha value that 
% produced the minimum sum of Zscore across our four fitting metrics.
%
% On further reflection it appears this is not quite the Z_score, which has
% no squared, but it is a measure of goodness of fit and it doesn't really
% change the results - it was an arbitrary choice of metric anyway, both do
% the job equally as well.
ZScore = zeros(length(alphas), 4);
for i = 1:4
    if i == 1
        X = squeeze(Accuracy(:,length(NumNeus),1));
        XE = squeeze(Deviation(:,length(NumNeus),1));
        Y = AccuracyData(length(NumNeus),1);
    elseif i == 2
        X = squeeze(AccuracySide(:,length(NumNeus)));
        XE = squeeze(DeviationSide(:,length(NumNeus)));
        Y = AccuracyDataSide(length(NumNeus));
    elseif i == 3
        X = CorrCoef;
        XE = CorrCoefDev;
        Y = CorrCoefData;
    else 
        X = BilCorr;
        XE = BilCorrDev;
        Y = BilCorrData;
    end
    for j = 1:length(alphas)
        ZScore(j, i) = ((X(j) - Y)^2)/XE(j);
    end
end
ZScore = sum(ZScore, 2);
[~, AlphaChoiceInd] = min(ZScore);
AlphaChoose = alphas(AlphaChoiceInd);

%% Now plotting, first comparison of ipsilateral macroscopic properties: Fig S7 A-C

figure
sgtitle('Plot 1: Ipsi Model and Data Comparison, Fig S7 A-C')

% Setup response magnitude distribution
hlim = -25:2:25;
subplot(1,3,1)
hold on
[Y, X1] = histcounts(abs(h1(h1~=0)),100,'Normalization','pdf');
BinWidth = X1(2) - X1(1);
X = X1(1:length(Y)) + BinWidth;
histogram(MagnitudesData{1},60,'Normalization','pdf', 'DisplayName', 'Data')
plot(X, Y, 'DisplayName', 'Model','linewidth',2)
legend
xlim([0,15])
xlabel('Response Magnitude (Hz)')
ylabel('PDF')
legend('Location','NorthWest')
hold off

% Next proportion responding to Odours
subplot(1,3,2)
hold on
Sums = sum(abs(h1)>0, 2);
histogram(OdoursData{2}, 'Normalization', 'probability', 'DisplayName', 'Data','BinEdges', OdourEdges)
[Y, X1] = histcounts(Sums,'Normalization','pdf','BinEdges', OdourEdges);
BinWidth = X1(2) - X1(1);
X = X1(1:length(Y)) + BinWidth/2;
plot(X, Y, 'DisplayName', 'Model','linewidth',2)
xlabel('Number of Odours')
ylabel('Proportion of Neurons')
legend
hold off

% Finally PNB calculations
% Here we make sure that points that fall exatly on a histogram bin
% boundary are shared equally between the two bins, rather than the matlab
% version that assigns them to the upper.
SumsP = sum(h1>0,2);
PNB = SumsP./Sums;
PNB2 = PNB(~isnan(PNB));
Box = zeros(5,1);
Total = length(PNB2);

for i = 1:5
    Box(i) = sum((PNB2 > 0.2*(i-1)) & (PNB2 < 0.2*i));
end

Box(1) = Box(1) + sum(PNB2 == 0);
Box(5) = Box(5) + sum(PNB2 == 1);

for i = 1:4
    Dividor = sum(PNB == 0.2*i);
    Box(i) = Box(i) + Dividor/2;
    Box(i+1) = Box(i+1) + Dividor/2;
end

Box = Box/Total;

% Before then plotting PNB
subplot(1,3,3)
hold on
histogram(PNBData{2},'Normalization', 'probability', 'DisplayName', 'Data','BinEdges', PNBEdges)
X = zeros(1,length(PNBEdges)-1);
for kappa = 1:(length(PNBEdges)-1)
    X(kappa) = PNBEdges(kappa)/2 + PNBEdges(kappa+1)/2;
end
plot(X, Box, 'DisplayName', 'Model','linewidth',2)
xlabel('Positive Negative Balance')
ylabel('Proportion of Responding Neurons')
legend('Location', 'North')
hold off

%% Next, how our four metrics of alignment comparison vary with alpha: Fig 6F,G,I and J
% Also Fig S7 E - H if you change the number of neurons simulated.
% Setup some useful things
xx = 0:0.01:1;
alphaplot = [1:8, AlphaChoiceInd, 9+9, 19:length(alphas)];
figure
sgtitle('Plot 2: Comparison of four alignment metrics vs alpha, Fig 6F,G,I and J')

% Percentage of Bilaterally Correlated Neurons
subplot(2,2,1)
z1 = [BilCorrData, BilCorrData];
z3 = [min(alphas),max(alphas)];
hold on
errorbar(alphas(alphaplot), BilCorr(alphaplot), BilCorrDev(alphaplot), '*', 'DisplayName', '385 Neurons')
errorbar(alphas(AlphaChoiceInd), BilCorr(AlphaChoiceInd), BilCorrDev(AlphaChoiceInd), 'r*', 'DisplayName', 'Best Alpha')
plot(z3, z1,'color',[0,0.5,0],'linewidth',0.5,'DisplayName', 'AON')
plot(alphas(alphaplot), BilCorr(alphaplot),'b','HandleVisibility','off')
title('Proportion of Bilaterally Correlated Neurons')
xlabel('\alpha')
ylabel('% Bil Corr')
legend('Location', 'SouthEast')
ylimits = [0, 0.7];
ylim(ylimits)
hold off

% Correlation Coefficient between Ipsi and Contra representations
subplot(2,2,2)
hold on
errorbar(alphas(alphaplot), CorrCoef(alphaplot), CorrCoefDev(alphaplot), '*','DisplayName', [num2str(385), ' Neurons'])
errorbar(alphas(AlphaChoiceInd), CorrCoef(AlphaChoiceInd), CorrCoefDev(AlphaChoiceInd), '*r','DisplayName','Best Alpha')
plot(alphas(alphaplot), CorrCoef(alphaplot),'b','HandleVisibility','off')
title('Correlation between all Ipsi and Contra Activities')
xlabel('\alpha')
ylabel('C vs I Correlation','DisplayName', 'Model Contra')
z1 = [CorrCoefData(1), CorrCoefData(1)];
z3 = [min(alphas),max(alphas)];
plot(z3, z1,'color',[0,0.5,0],'linewidth',0.5,'DisplayName', 'AON')
legend('Location','SouthEast')
ylimits = [-0.01, 0.8];
ylim(ylimits)
hold off


Index = 13;
% Decoding work, first decoding odours
subplot(2,2,3)
hold on
errorbar(alphas(alphaplot), Accuracy((alphaplot),Index,1),Deviation((alphaplot),Index,1), '*', 'DisplayName', [num2str(385), ' Neurons'])
errorbar(alphas(AlphaChoiceInd), Accuracy((AlphaChoiceInd),Index,1),Deviation((AlphaChoiceInd),Index,1), '*r', 'DisplayName','Best Alpha')
plot(alphas(alphaplot), Accuracy((alphaplot),Index,1),'b','HandleVisibility','off')
linex = min(alphas):0.01:max(alphas);
lineyup = repmat(max(AccuracyData(Index,1)), [1, length(linex)]);
legend('Location', 'SouthEast')
plot(linex, lineyup,'color',[0,0.75,0],'linewidth',0.5,'DisplayName', 'AON')
xlabel('\alpha')
ylabel('Accuracy')
ylimits = [0,1];
ylim(ylimits)
title('Accuracy - Classifying Contralateral Odor Centroid IDs using Ipsilateral Centroids')
hold off

% Finally decoding sides
subplot(2,2,4)
hold on
Index = length(NumNeus);
errorbar(alphas(alphaplot), AccuracySide((alphaplot),Index),DeviationSide((alphaplot),Index),'*', 'DisplayName', [num2str(385), ' Neurons'])
errorbar(alphas(AlphaChoiceInd), AccuracySide((AlphaChoiceInd),Index),DeviationSide((AlphaChoiceInd),Index),'*r', 'DisplayName', 'Best Alpha')
linex = min(alphas):0.01:max(alphas);
lineyup = repmat(max(AccuracySide(:,Index)), [1, length(linex)]);
plot(linex, lineyup,'color',[0,0.75,0],'linewidth',0.5,'DisplayName', 'AON')
plot(alphas(alphaplot), AccuracySide((alphaplot),Index),'b','HandleVisibility','off')
xlabel('alpha')
ylabel('Accuracy')
title('Accuracy - Centroid Side Classification Task')
ylimits = [0,1];
ylim(ylimits)
legend('Location','South')
hold off

%% How well contra reps can decode ipsi depends on the decoder: fig S7D
figure
sgtitle('Plot 3: Contra Decoding Ipsi Discussion Fig S7D')

% First up we plot the simple euclidean decoder
subplot(1,2,1)
hold on
errorbar(alphas(alphaplot), Accuracy((alphaplot),Index,2),Deviation((alphaplot),Index,2), '*', 'DisplayName', [num2str(385), ' Neurons'])
errorbar(alphas(AlphaChoiceInd), Accuracy((AlphaChoiceInd),Index,2),Deviation((AlphaChoiceInd),Index,2), '*r', 'DisplayName', 'Best Alpha')
linex = min(alphas):0.01:max(alphas);
lineyup = repmat(max(AccuracyData(Index,2)), [1, length(linex)]);
yy = spline(alphas, Accuracy(:,Index,2), xx);
plot(xx, yy,'b','HandleVisibility','off')
legend('Location', 'SouthEast')
plot(linex, lineyup,'color',[0,0.75,0],'linewidth',0.5,'DisplayName', 'AON')
xlabel('alpha')
ylabel('Proportion of Centroids Correctly Classified')
title('Contra Decoding Ipsi - Euclidean Decoding')
legend
hold off

% Then we show the angular one, hence showing it is length discrepency not
% lack of information that is stalling the euclidean decoder.
subplot(1,2,2)
hold on
errorbar(alphas(alphaplot), AccuracyAng(alphaplot),DeviationAng(alphaplot), '*', 'DisplayName', [num2str(385), ' Neurons'])
errorbar(alphas(AlphaChoiceInd), AccuracyAng(AlphaChoiceInd),DeviationAng(AlphaChoiceInd), '*r', 'DisplayName', 'Best Alpha')
linex = min(alphas):0.01:max(alphas);
lineyup = repmat(AccuracyData2, [1, length(linex)]);
yy = spline(alphas, AccuracyAng, xx);
plot(xx, yy,'b','HandleVisibility','off')
legend('Location', 'SouthEast')
points = [0.36, 0.46];
ylim([0,1])
plot(linex, lineyup,'color',[0,0.75,0],'linewidth',0.5,'DisplayName', 'AON')
ylim([0,1])
xlabel('alpha')
ylabel('Proportion of Centroids Correctly Classified')
title('Contra Decoding Ipsi - Angle Decoding')
legend
hold off

%% Now we compare the contralateral macro distributions to those in data, Fig 6B-D
figure
sgtitle('Plot 4: Comparison of Contralateral simulation and data representation, Fig 6B-D')

alphachoose = 15;
offset = 0;
Comps = length(alphachoose)+ offset;
plots = 3;
maglims = [0,15];
maglims2 = [-12, 12];
maglimy = [0, 0.22];

% Plot the alpha varying things
for  i = 1:length(alphachoose)
    param = alphas(alphachoose(i));
    resp = z(:,:,alphachoose(i));
    
    subplot(Comps, plots,1+plots*(i-1))
    hold on
    [Y, X1] = histcounts(abs(resp(resp~=0)),60,'Normalization','pdf');
    BinWidth = X1(2) - X1(1);
    X = X1(1:length(Y)) + BinWidth/2;
    histogram(MagnitudesData{1},60,'Normalization','pdf', 'DisplayName', 'Data')
    plot(X, Y, 'DisplayName', 'Model','linewidth',2)
    legend
    xlim([0,15])
    xlabel('Response Magnitude (Hz)')
    ylabel('PDF')
    
    subplot(Comps, plots, 2+plots*(i-1))
    hold on
    Sums = sum(abs(resp)>0, 2);
    histogram(OdoursData{1}, 'Normalization', 'probability', 'DisplayName', 'Data','BinEdges', OdourEdges)
    [Y, X1] = histcounts(Sums,'Normalization','pdf','BinEdges',OdourEdges);
    BinWidth = X1(2) - X1(1);
    X = X1(1:length(Y)) + BinWidth/2;
    plot(X, Y, 'DisplayName', 'Model','linewidth',2)
    xlabel('Odours')
    ylabel('Proportion of Neurons')
    hold off
    
    subplot(Comps, plots, 3+plots*(i-1))
    % Plot the bestest PNB
    SumsP = sum(resp>0,2);
    PNB = SumsP./Sums;
    PNB2 = PNB(~isnan(PNB));
    Box = zeros(5,1);
    Total = length(PNB2);

    for k = 1:5
        Box(k) = sum((PNB2 > 0.2*(k-1)) & (PNB2 < 0.2*k));
    end

    Box(1) = Box(1) + sum(PNB2 == 0);
    Box(5) = Box(5) + sum(PNB2 == 1);

    for k = 1:4
        Dividor = sum(PNB == 0.2*k);
        Box(k) = Box(k) + Dividor/2;
        Box(k+1) = Box(k+1) + Dividor/2;
    end

    Box = Box/Total;

    subplot(1,3,3)
    hold on
    histogram(PNBData{2},'Normalization', 'probability', 'DisplayName', 'Data','BinEdges', PNBEdges)
    X = zeros(1,length(PNBEdges)-1);
    for kappa = 1:(length(PNBEdges)-1)
        X(kappa) = PNBEdges(kappa)/2 + PNBEdges(kappa+1)/2;
    end
    plot(X, Box, 'DisplayName', 'Model','linewidth',2)
    xlabel('Positive Negative Balance')
    ylabel('Proportion of Responding Neurons')
    legend('Location', 'North')
    hold off
end

%% Here we compare a random and optimally structured ipsi vs contra data cloud, Fig 6E
figure
alphachoose = [1, AlphaChoiceInd];
sgtitle('Plot 5: Comparison of random and optimally structured ipsi vs contra data cloud Fig 6E')

for i = 1:length(alphachoose)
    param = alphas(alphachoose(i));
    resp = z(:,:,alphachoose(i));
    subplot(1, length(alphachoose), i)
    hold on
    xlabel('Ispi Response (Hz)')
    ylabel('Contra Response (Hz)')
    Corker = corr(h1(:), resp(:));
    p1 = scatter(h1(:), resp(:), 'filled');
    p1.MarkerFaceAlpha = 0.1;
    b1 = h1(:)\resp(:);
    yCalc = b1*h1(:);
    plot(h1(:), yCalc)
    xlim(maglims2)
    ylim(maglims2)
    txt = sprintf(['Grad: ', num2str(round(b1,2,'significant')),'\n','Corr: ', num2str(round(Corker,2,'significant'))]);
    text(2,-7,txt)
end
%% Graph showing decoding ability with neural population size, random vs optimally structured, Fig 6H

figure
title('Plot 6: Decoding ability with population size, Fig 6H')

hold on
xlabel('Number of Neurons')
ylabel('Accuracy')
alphachoose = [1,AlphaChoiceInd];
for i = 1:length(alphachoose)
    errorbar(NumNeus, Accuracy(alphachoose(i),:,1), Deviation(alphachoose(i),:,1),'DisplayName', ['\alpha: ', num2str(alphas(alphachoose(i)))])
end
errorbar(NumNeus, AccuracyData(:,1), DeviationData(:,1), 'DisplayName', 'Data')
z1 = [0, NumNeuExp];
z2 = [1/15, 1/15];
plot(z1, z2, 'k', 'DisplayName', 'Chance')
legend('Location', 'East')
hold off