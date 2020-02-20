%% Data Extractor.m
% Our later modelling work will be validate through comparison to data.
% Here we extract from the original data the required parameters and plots
% for comparison. These will then be saved in
% NeuronalDataforModelComparison.mat
%
% Since the model only concerns AON we do not use APC or PPC data.
% However the model takes the cosine similarity of the OB odour
% respresentations as input, therefore we also extract that.

%% First up we load the data and extract the 
load('Data_and_Notes\tetrodeRecordings_OC_2s.mat')
clear BasalFR B C
alpha = 0.05;
Neurons = 385;
Odours = 15;
Trials = 7;
Signif = zeros(Neurons,Odours,2);
SignifError = zeros(Neurons, Odours, 2, Trials);
Sides = {'Contra', 'Ipsi'};

AON = cat(1,A{1:3});
Base = squeeze(AON(:,16,:,:));
AON = AON(:,1:Odours,:,:);

% Here we calculate whether the deviation of the data from basal rate was
% significant (by a signed rank test)
for side = 1:2
    BaseAvg = mean(squeeze(Base(:,side,:)),2);
    for n = 1:Neurons
        for o = 1:Odours
            Signif(n, o, side) = signrank(squeeze(AON(n, o, side, :)), BaseAvg(n))<=alpha;
            for t = 1:Trials
                List = 1:7;
                List(t) = [];
                SignifError(n, o, side, t) = signrank(squeeze(AON(n, o, side, List)), BaseAvg(n))<=alpha;
            end
        end
    end
    % Before finally subtracting the base
    AON(:,:,side,:) = AON(:,:,side,:) - repmat(BaseAvg,[1,Odours,1,Trials]);
end

% We only compare to the trial averaged data
AON_Avg = mean(AON, 4);

%% Here we plot and extract from the AON distributions to which we shall compare
% This code plots a version of the top (AON) graphs from Fig 3E-G.
%
% NOT PARTICULARLY IMPORTANT COMMENT
% These will not look identical to the ones in the paper (despite the fact
% that you rightly think they should!). This is due to the way matlab plots
% histograms. Say there is one bin covering the range [0.6, 0.8] and another
% [0.8, 1], where should you put 0.8?
%
% Matlab chooses the upper bin. This can make a fairly symmetric
% distribution look assymmetric. In the paper we assign half to the upper
% bin and half to the lower. Since these plots created here are just to
% make sure it looks okay before we then extract all the data for later
% comparison we don't worry about it and just extract the data.
%
% Odour_X_ correponds to fig 3E, Mag_X_ to 3F, PNB_X_ to 3G.
OdourEdges = 0:1:15;
PNBEdges = 0:0.2:1;
MagEdges = 0:2:25;
       
MagnitudesData = cell(2,1);
PNBData = MagnitudesData;
PNBPlot = MagnitudesData;
OdoursData = MagnitudesData;
OdourPlot = MagnitudesData;
OdourDataError = zeros(Neurons, Trials, 2);
PNBDataError = MagnitudesData;
OdoursVar = zeros(length(OdourEdges)-1, 2);
PNBVar = zeros(length(PNBEdges) - 1, 2);
S = zeros(2, 1);
MagX = MagnitudesData;
MagY = MagnitudesData;

% Start looking at key plots
figure
for side = 1:2
    AON_Side = AON_Avg(:,:,side);
    Signif_Side = Signif(:,:,side);
    SignifError_Side = squeeze(SignifError(:,:,side,:));
    
    % Extract Magnitudes
    MagnitudesData{side} = abs(AON_Side(Signif_Side(:,:)~=0));
    S(side) = sum(Signif_Side, 'all')/numel(Signif_Side);
    subplot(2,3,1+(side-1)*3)
    Mag = histogram(MagnitudesData{side},MagEdges, 'Normalization', 'pdf','DisplayName',[Sides{side},' Data']);
    MagY{side} = Mag.Values;
    Mag_Here = MagY{side};
    for i = 1:length(Mag_Here)
        Mag_Here(i) = (MagEdges(i+1)+MagEdges(i))/2;
    end
    MagX{side} = Mag_Here;
    title(['Distribution of Response Magnitudes - ', Sides{side}])
    ylabel('Probability')
    xlabel('Response (Hz)')

    % Now odours
    subplot(2,3,2+(side-1)*3)
    OdoursData{side} = sum(Signif_Side~=0, 2);
    OdoursVartmp = zeros(length(OdourEdges)-1, Trials);
    for t = 1:Trials
        OdourDataError(:,t,side) = sum(SignifError_Side(:,:,t)~=0, 2);
        Hist =  histogram(OdourDataError(:,t,side), 'Normalization', 'probability','BinEdges', OdourEdges);
        OdoursVartmp(:,t) = Hist.Values;
    end
    OdoursVar(:,side) = std(OdoursVartmp, 0, 2);
    Hist = histogram(OdoursData{side}, 'Normalization', 'probability','BinEdges', OdourEdges,'DisplayName',[Sides{side},' Data']);
    hold on
    errorbar(OdourEdges(1:Odours)+0.5, Hist.Values, OdoursVar(:,side), 'k', 'linestyle', 'none','DisplayName','6 Trial Error')
    OdourPlot{side} = Hist.Values;
    title(['Distribution of Odour Responders - ', Sides{side}])
    ylabel('Proportion of Neurons')
    xlabel('Number of Odours')
    
    % Now PNB
    subplot(2,3,3 + (side-1)*3)
    OdoursData_Side = OdoursData{side};
    Positive = sum((Signif_Side~=0 & AON_Side > 0), 2);
    PNBData{side} = Positive(OdoursData_Side~=0)./OdoursData_Side(OdoursData_Side~=0);
    PNBDataError_Side = cell(Trials, 2);
    PNBVartmp = zeros(length(PNBEdges)-1, Trials);
    for t = 1:Trials
        Positive = sum((SignifError_Side(:,:,t)~=0 & AON_Side > 0), 2);
        OdourDataError_Here = OdourDataError(:,t,side);
        PNBDataError_Side{t,side} = Positive(OdourDataError_Here~=0)./OdourDataError_Here(OdourDataError_Here~=0);
        Hist =  histogram(PNBDataError_Side{t,side}, 'Normalization', 'probability','BinEdges', PNBEdges);
        PNBVartmp(:, t) = Hist.Values;
    end
    PNBVar(:,side) = std(PNBVartmp, 0, 2);
    PNBDataError{side} = PNBDataError_Side;
    Hist = histogram(PNBData{side}, 'Normalization', 'probability','BinEdges', PNBEdges,'DisplayName',[Sides{side},' Data']);
    hold on
    errorbar(PNBEdges(1:length(PNBEdges)-1)+0.1, Hist.Values, PNBVar(:,side), 'k', 'linestyle', 'none','DisplayName','6 Trial Error')
    PNBPlot{side} = Hist.Values;
    title(['Positive Negative Balance - ', Sides{side}])
    ylabel('Proportion of Responding Neurons')
    xlabel('PNB')
end

sgtitle('Plot 1: Checking to see if distributions look alright')
clear A BaseAvg cutoff i Mag Negative Positive side Signif_Side AON_Side Sums Hist List Mag_Here n o OdourDataError_Here OdoursData_Side OdoursVartmp PNBDataError_Side PNBVartmp SignifError_Side t

%% Now fit the sparsenss and magnitude scaling factor for Ipsi and Contra]
% Next we try and fit a tail of gaussian distribution to the response
% magnitudes. We use the best fitting scale conversion as a way to convert
% all numbers to firing rates later on.
MagFitI = @(y)ErrorMag(y,MagX{2},MagY{2},S(2));
x0 = 1;
ScaleConversion = fminsearch(MagFitI,x0);

% Plotting Section to check it looks okay
phioversigma = sqrt(2)*erfcinv(S(2));
Fit = @(x) sqrt(1/(2*pi))*ScaleConversion/S(2).*exp(-0.5.*(ScaleConversion.*x+sign(x).*phioversigma).^2);
subplot(2,3,4)
resolution = 0.1;
X = [0:resolution:-resolution,resolution:resolution:15];
yfit = Fit(X);
bar(MagX{2},MagY{2}, 1,'DisplayName',[Sides{2},' Data']);
hold on
plot(X,yfit,'r','DisplayName', 'Fitting');
xlabel('Magnitude of Response (Hz)')
ylabel('Probability')
title(['Distribution of Response Magnitudes - ', Sides{2}])
legend

%% Finally we load the OB data in order to extract the correlation matrix.
load('Data_and_Notes\tetrode_recordings_formatted_for_alex_v2_ordered_ob.mat')
load('Data_and_Notes\SniffTimeOB.mat')

ArrayOB = [ArrayOB1; ArrayOB2; ArrayOB3];
[OBData, ~] = DataCreator(1, 3, 1:15, 2, 1, {ArrayOB}, SniffTimeOB);

OBMean = mean(OBData, 3);
CosineSim = squareform(getCosineSimilarity(OBMean)) + eye(15);
figure
imagesc(CosineSim)
%%
save('Data_and_Notes\NeuronalDataforModelComparison.mat','OdoursData','PNBData','MagnitudesData','S','ScaleConversion', 'PNBVar', 'OdoursVar', ...
    'PNBEdges', 'MagY', 'MagX', 'OdourEdges', 'OdourPlot', 'PNBPlot', 'MagEdges','AON', 'CosineSim')
