function code_for_figure_2()

%%%%%%%%%%%%%%%%%%%
%FIGURES 2C AND 2D%
%%%%%%%%%%%%%%%%%%%

load('gcamp3_mouse1_intensityTraces_glomerularMasks')

%Parameters
nbOdors = 16;
odorToDisplay = 1;
glomToDisplay = [1 2];
fps = 4; %video captured at 4 frames per seconde
baseline = 4; %in s
nbTrials = 7;

%Calculate delta f over f for each odor, each side, each glomerulus
s = size(intensityTracesL);
deltaFoverF = zeros(2,nbOdors,s(3));
for G = 1:s(3)
    for O = 1:nbOdors
        basalFL = mean(mean(intensityTracesL(O:16:16*7,1:baseline*fps,G)));
        responseFL = max(mean(intensityTracesL(O:16:16*7,baseline*fps+1:s(2),G)));
        deltaFoverF(1,O,G) = (responseFL - basalFL) / basalFL;

        basalFR = mean(mean(intensityTracesR(O:16:16*7,1:baseline*fps,G)));
        responseFR = max(mean(intensityTracesR(O:16:16*7,baseline*fps+1:s(2),G)));
        deltaFoverF(2,O,G) = (responseFR - basalFR) / basalFR;
    end
end

%Compute heat maps of the bulb
GlomeruliInvert = 1 - Glomeruli;
imageL = zeros(1944,2592,nbOdors);
imageR = zeros(1944,2592,nbOdors);
for G = 1:s(3)
    for O = 1:nbOdors
        imageL(:,:,O) = imageL(:,:,O) + GlomeruliInvert(:,:,G) * deltaFoverF(1,O,G);
        imageR(:,:,O) = imageR(:,:,O) + GlomeruliInvert(:,:,G) * deltaFoverF(2,O,G);
    end
end

%Display heat map + picture of the craniotomy (Figure 2C)
Mask = zeros(1944,2592);
for g = 1:s(3)
    Mask = Mask + edge(Glomeruli(:,:,g));
end
Mask(Mask>1) = 1;
Mask = imdilate(Mask, strel('disk',1));
Mask = imdilate(Mask, strel('disk',1));
figure
subplot(3,1,1)
IM = mipbaseline.*(-1*(Mask-1))+Mask*211;
imagesc(IM(575:1400,475:1950),[0 180])
c = hot(100);
mini = 0;
maxi = 12/100;
subplot(3,1,2)
h = imagesc(imageL(575:1400,475:1950,odorToDisplay),[mini maxi]);
set(h,'AlphaData',imageL(575:1400,475:1950,odorToDisplay)~=0);
colormap(c(1:70,:))
subplot(3,1,3)
h = imagesc(imageR(575:1400,475:1950,odorToDisplay),[mini maxi]);
set(h,'AlphaData',imageR(575:1400,475:1950,odorToDisplay)~=0);
colormap(c(1:70,:))

%Display examplar glomeruli (Figure 2D)
for G = glomToDisplay
    LL = zeros(7,37);
    for t = 1:nbTrials
        LL(t,:) = (intensityTracesL(odorToDisplay +(t-1)*16,:,G) - mean(intensityTracesL(odorToDisplay +(t-1)*16,1:baseline*fps,G),2))/mean(intensityTracesL(odorToDisplay +(t-1)*16,1:baseline*fps,G),2);
    end
    RR = zeros(7,37);
    for t = 1:nbTrials
        RR(t,:) = (intensityTracesR(odorToDisplay +(t-1)*16,:,G) - mean(intensityTracesR(odorToDisplay +(t-1)*16,1:baseline*fps,G),2))/mean(intensityTracesR(odorToDisplay +(t-1)*16,1:baseline*fps,G),2);
    end
    figure
    hold on
    LLmean = mean(LL);
    RRmean = mean(RR);
    LLsem = std(LL)/sqrt(7);
    RRsem = std(RR)/sqrt(7);
    shadedErrorBar([-4:1/fps:5],LLmean,LLsem,'lineprops',{'b','markerfacecolor','b'})
    shadedErrorBar([-4:1/fps:5],RRmean,RRsem,'lineprops',{'r','markerfacecolor','r'})
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%
%FIGURES 2E AND 2H%
%%%%%%%%%%%%%%%%%%%

load('gcamp3_allmice_glomResponsesIpsiSignif')
figure

%Plot the contra vs ipsi responses (Figures 2E and 2H)
X = [dFFOneSide_m1_b1';dFFOneSide_m1_b2';dFFOneSide_m2_b1';dFFOneSide_m2_b2';dFFOneSide_m3_b1';dFFOneSide_m3_b2'];
Y = [dFFOneSide2_m1_b1';dFFOneSide2_m1_b2';dFFOneSide2_m2_b1';dFFOneSide2_m2_b2';dFFOneSide2_m3_b1';dFFOneSide2_m3_b2'];
subplot(1,2,1)
hold on
plot(X,Y,'ok')
plot([-0.05 0.25],[-0.05 0.25],'k')
plot(mean(X),mean(Y),'or')
xlim([-0.02 0.25])
ylim([-0.02 0.25])
grid on
axis square
[mean(X) mean(Y) mean(Y)/mean(X)]

%Plot the histogram for contra responses, fit a Gaussian curve (Figure 2E)
c = 0.005;
stp = [-0.02:c:0.25];
subplot(1,2,2)
hold on
h = histogram(Y,stp);
stpG = [-0.02+c/2:c:0.25];
h = h.Values;
[f,gof] = fit(stpG.',h.','gauss1');
plot(f,stpG,h)
xlim([-0.02 0.25])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%
%FIGURE 2G%
%%%%%%%%%%%

load('tetrodeRecordings_OB_new')

%Plot nb of units responding to ipsi and/or contra stimulations (Figure 2G)
M = C{1};
M = M(:,1:15,:); %the 16th odor is the blank trial
Mipsi = sum(M(:,:,2),2)&1;
Mcontra = sum(M(:,:,1),2)&1;
Ionly = sum(Mipsi&~Mcontra);
Conly = sum(Mcontra&~Mipsi);
IandC = sum(Mipsi&Mcontra);
noR = size(M,1) - Ionly - Conly - IandC;
figure
plot([Ionly Conly IandC noR])

%Plot pie chart for percentage positive / negative responses (Figure 2G)
M = C{1};
M = M(:,1:15,:); %the 16th odor is the blank trial
R = A{1};
R = R(:,1:15,:);
Ripsi = R(:,:,2).*M(:,:,2);
Rcontra = R(:,:,1).*M(:,:,1);
howManyPos_Ripsi = length(find(Ripsi>0));
howManyNeg_Ripsi = length(find(Ripsi<0));
howManyPos_Rcontra = length(find(Rcontra>0));
howManyNeg_Rcontra = length(find(Rcontra<0));
figure
subplot(1,2,1)
pie([howManyPos_Ripsi howManyNeg_Ripsi])
subplot(1,2,2)
pie([howManyPos_Rcontra howManyNeg_Rcontra])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%
%FIGURE 2I%
%%%%%%%%%%%

%Plot the contra vs ipsi responses (Figures 2I)
load('tetrodeRecordings_OB_new')
M = C{1};
M = M(:,1:15,:); %the 16th odor is the blank trial
R = A{1};
R = abs(R(:,1:15,:));
Mipsi = M(:,:,2);
Ripsi = R(:,:,2);
Rcontra = R(:,:,1);
Ripsi(find(Mipsi==0)) = NaN;
Rcontra(find(Mipsi==0)) = NaN;
figure
hold on
plot(Ripsi,Rcontra,'ok')
plot(nanmean(reshape(Ripsi,42*15,1)),nanmean(reshape(Rcontra,42*15,1)),'or')
xlim([0 6])
ylim([0 6])


