function code_for_figure_4()

%%%%%%%%%%%
%FIGURE 4B%
%%%%%%%%%%%

Percentages = zeros(10,1);
alpha = 0.05;

%Percentage of bilaterally-correlated neurons in awake mice (Figure 4B)
load('tetrodeRecordings_OC_2s.mat')
for mouse = 1:10
    mouse
    PVAL = [];
    M = A{mouse};
    M = mean(M,4);
    for odor = 1:15
        M(:,odor,:) = M(:,odor,:) - M(:,16,:);
    end
    Mi = M(:,1:15,2);
    Mc = M(:,1:15,1); 
    for neuron = 1:size(Mi,1)
        MNi = Mi(neuron,:);
        MNc = Mc(neuron,:);
        lm = fitlm(MNi,MNc,'linear');
        lm = lm.Coefficients;
        lm = table2array(lm);
        pval = lm(2,4);
        PVAL = [PVAL;pval];
    end
    Percentages(mouse) = length(find(PVAL<=alpha))*100/length(PVAL);
    [length(find(PVAL<=alpha)) length(PVAL) length(find(PVAL<=alpha))*100/length(PVAL)]
end

%Display percentage of bilaterally-correlated neurons per mouse (Figure 4B)
figure
plot([1 1 1 2 2 2 2 3 3 3],Percentages,'o')
ylim([0 50])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%
%FIGURE 4C%
%%%%%%%%%%%

%Percentage of bilaterally-corr neurons, bootstrap analysis (Figure 4C)
alpha = 0.05;
load('tetrodeRecordings_OC_2s.mat')
distr = [];
mouse = 1; %repeat for all mice
nbRepBootstrap = 1000;
for i = 1:nbRepBootstrap
    i*100/nbRepBootstrap
    PVAL = [];
    M = A{mouse};
    M = mean(M,4);
    for odor = 1:15
        M(:,odor,:) = M(:,odor,:) - M(:,16,:);
    end
    Mi = M(:,1:15,2);
    Mc = M(:,1:15,1);
    for neuron = 1:size(Mi,1)
        MNi = Mi(neuron,:);
        MNc = Mc(neuron,:);
        MNi = MNi(randperm(15));
        MNc = MNc(randperm(15));
        lm = fitlm(MNi,MNc,'linear');
        lm = lm.Coefficients;
        lm = table2array(lm);
        pval = lm(2,4);
        PVAL = [PVAL;pval];
    end
    percentage = length(find(PVAL<=alpha))*100/length(PVAL);
    distr = [distr;percentage];
end
figure
histogram(distr,0:2:max(distr))
pvalForBootstrap = length(find(distr>Percentages(mouse)))/length(distr);
title(['Bootstrap mean:' num2str(mean(distr)) '. Actual mean:' num2str(Percentages(mouse)) '. P-value:' num2str(pvalForBootstrap)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%
%FIGURES 4D, 4E, AND 4F%
%%%%%%%%%%%%%%%%%%%%%%%%

%R2, slope, and intercept of bilat-corr neurons in the AON (Figures 4D-F)
RSI = [];
alpha = 0.05;
load('tetrodeRecordings_OC_2s.mat')
for mouse = 1:7; %includes all AON and APC mice from awake recordings
    mouse
    PVAL = [];
    M = A{mouse};
    M = mean(M,4);
    for odor = 1:15
        M(:,odor,:) = M(:,odor,:) - M(:,16,:);
    end
    Mi = M(:,1:15,2);
    Mc = M(:,1:15,1);
    for neuron = 1:size(Mi,1)
        MNi = Mi(neuron,:);
        MNc = Mc(neuron,:);
        lm = fitlm(MNi,MNc,'linear');
        lm = lm.Coefficients;
        lm = table2array(lm);
        pval = lm(2,4);
        if pval <= alpha
            p = polyfit(MNi,MNc,1);
            yfit =  p(1) * MNi + p(2);
            yresid = MNc - yfit;
            SSresid = sum(yresid.^2);
            SStotal = (length(MNc)-1) * var(MNc);
            rsq = 1 - SSresid/SStotal;
            RSI = [RSI;rsq p(1) p(2)];
        end
    end
end
figure
Hlim = {};
Hlim{1} = [0:0.1:1];
Hlim{2} = [-1:0.2:1];
Hlim{3} = [-4:0.4:4];
for i = 1:3
    subplot(3,1,i)
    histogram(RSI(:,i),Hlim{i},'DisplayStyle','stairs')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%
%FIGURE 4A%
%%%%%%%%%%%

load('tetrodeRecordings_OC_2s.mat')

%Examplar bilaterally-correlated neuron from AON1 (Figure 4A)
mouse = 1;
M = A{mouse};
M = mean(M,4);
for odor = 1:15
    M(:,odor,:) = M(:,odor,:) - M(:,16,:);
end
Mi = M(:,1:15,2);
Mc = M(:,1:15,1);
count = 1;
for neuron = [3 17 2 83   10 114 185 47]
    MNi = Mi(neuron,:);
    MNc = Mc(neuron,:);
    lm = fitlm(MNi,MNc,'linear');
    lm = lm.Coefficients;
    lm = table2array(lm);
    pval = lm(2,4)
    p = polyfit(MNi,MNc,1);
    yfit =  p(1) * MNi + p(2);
    yresid = MNc - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(MNc)-1) * var(MNc);
    rsq = 1 - SSresid/SStotal;
    RSI = [rsq p(1) p(2)];
    
    subplot(2,4,count)
    hold on
    plot(MNi,MNc,'.k')
    plot([-5 5],[-5 5]*p(1) + p(2),'-r')
    axis square
    grid on
    xlim([-10 10])
    ylim([-10 10])
    title(RSI)
    count = count+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%
%FIGURE 3H%
%%%%%%%%%%%

%Odor decoding per region (Figure 3H)
cd('machine learning results')
load('ByCortexData31July.mat')
cd ../
figure
X = [0:10:50 100:50:384 384];
Y = Accuracies{1};
E = Deviations{1};
subplot(3,2,1)
hold on
shadedErrorBar(X,Y(4,:),E(4,:),'lineprops',{'r-o','markerfacecolor','k','markersize',2})
shadedErrorBar(X,Y(3,:),E(3,:),'lineprops',{'r-.x','markerfacecolor','k'})
plot([0 384],[1 1]/15,'-g')
subplot(3,2,2)
hold on
shadedErrorBar(X,Y(1,:),E(1,:),'lineprops',{'r-o','markerfacecolor','k','markersize',2})
shadedErrorBar(X,Y(2,:),E(2,:),'lineprops',{'r-.x','markerfacecolor','k'})
plot([0 384],[1 1]/15,'-g')

X = [0:10:50 100:50:931 931];
Y = Accuracies{2};
E = Deviations{2};
subplot(3,2,3)
hold on
shadedErrorBar(X,Y(4,:),E(4,:),'lineprops',{'b-o','markerfacecolor','k','markersize',2})
shadedErrorBar(X,Y(3,:),E(3,:),'lineprops',{'b-.x','markerfacecolor','k'})
plot([0 931],[1 1]/15,'-g')
subplot(3,2,4)
hold on
shadedErrorBar(X,Y(1,:),E(1,:),'lineprops',{'b-o','markerfacecolor','k','markersize',2})
shadedErrorBar(X,Y(2,:),E(2,:),'lineprops',{'b-.x','markerfacecolor','k'})
plot([0 931],[1 1]/15,'-g')

X = [0:10:50 100:50:505 505];
Y = Accuracies{3};
E = Deviations{3};
subplot(3,2,5)
hold on
shadedErrorBar(X,Y(4,:),E(4,:),'lineprops',{'k-o','markerfacecolor','k','markersize',2})
shadedErrorBar(X,Y(3,:),E(3,:),'lineprops',{'k-.x','markerfacecolor','k'})
plot([0 505],[1 1]/15,'-g')
subplot(3,2,6)
hold on
shadedErrorBar(X,Y(1,:),E(1,:),'lineprops',{'k-o','markerfacecolor','k','markersize',2})
shadedErrorBar(X,Y(2,:),E(2,:),'lineprops',{'k-.x','markerfacecolor','k'})
plot([0 505],[1 1]/15,'-g')

titles = {'AON train=ipsi','AON train=contra','APC train=ipsi','APC train=contra','PPC train=ipsi','PPC train=contra'};
for i = 1:6
    subplot(3,2,i)
    ylim([0 1])
    xlabel('number of neurons')
    ylabel('accuracy')
    title(titles{i});
    grid on
end

for region = 1:3
    for condition = 1:4
        X = squeeze(AccSpread(region,condition,:));
        signrank(X,1/15)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%
%FIGURES 4G & S4%
%%%%%%%%%%%%%%%%

%Ipsi versus contra responses per mouse, per odor (Figures 4G & S4)
close all
load('tetrodeRecordings_OC_2s.mat')
whichMouse = [1 1 1 2 2 2 2 3 3 3];
VALS = zeros(3,15,4);
for region = 1:3
    figure
    for odor = 1:15
        X = [];
        Y = [];
        for mouse_id = find(whichMouse==region)
            Mouse = A{mouse_id};
            Blank = mean(squeeze(Mouse(:,16,:,:)),3);
            Rall = mean(squeeze(Mouse(:,odor,2,:)),2);
            Lall = mean(squeeze(Mouse(:,odor,1,:)),2);
            Rall = Rall - Blank(:,2);
            Lall = Lall - Blank(:,1);
            X = [X;Rall];
            Y = [Y;Lall];
        end
%         X = X(randperm(length(X)));
        lm = fitlm(X,Y,'linear');
        lm = lm.Coefficients;
        lm = table2array(lm);
        pval = lm(2,4);
        
        Rall = X;
        Lall = Y;
        subplot(4,4,odor)
        hold on
        plot(Rall,Lall,'.k')
        p = polyfit(Rall,Lall,1);
        yfit =  p(1) * Rall + p(2);
        yresid = Lall - yfit;
        SSresid = sum(yresid.^2);
        SStotal = (length(Lall)-1) * var(Lall);
        rsq = 1 - SSresid/SStotal;
        VALS(region,odor,:) = [pval p(1) p(2) rsq];
        lm = 15;
        plot([-lm lm],[-lm lm]*p(1)+p(2),'r')
        xlim([-lm lm])
        ylim([-lm lm])
        grid on
        axis square
    end
end

val = 1;
M = VALS(:,:,val)';
y = repmat(1:15,3,1)'; % generate x-coordinates
x = repmat(1:3,15,1); % generate y-coordinates
% Generate Labels
t = num2cell(M); % extact values into cells
t = cellfun(@num2str, t, 'UniformOutput', false); % convert to string
% Draw Image and Label Pixels
figure
imagesc(M)
text(x(:), y(:), t, 'HorizontalAlignment', 'Center')


%Ipsi versus contra responses per mouse, per odor (Figure S4)
load('tetrodeRecordings_OC_2s.mat')
VALS = zeros(10,15,4);
for mouse_id = 1:10
    figure
    for odor = 1:15
        Mouse = A{mouse_id};
        Blank = mean(squeeze(Mouse(:,16,:,:)),3);
        Rall = mean(squeeze(Mouse(:,odor,2,:)),2);
        Lall = mean(squeeze(Mouse(:,odor,1,:)),2);
        X = Rall - Blank(:,2);
        Y = Lall - Blank(:,1);
        lm = fitlm(X,Y,'linear');
        lm = lm.Coefficients;
        lm = table2array(lm);
        pval = lm(2,4);
        
        Rall = X;
        Lall = Y;
        subplot(4,4,odor)
        hold on
        plot(Rall,Lall,'.k')
        p = polyfit(Rall,Lall,1);
        yfit =  p(1) * Rall + p(2);
        yresid = Lall - yfit;
        SSresid = sum(yresid.^2);
        SStotal = (length(Lall)-1) * var(Lall);
        rsq = 1 - SSresid/SStotal;
        VALS(mouse_id,odor,:) = [pval p(1) p(2) rsq];
        lm = 15;
        plot([-lm lm],[-lm lm]*p(1)+p(2),'r')
        xlim([-lm lm])
        ylim([-lm lm])
        grid on
        axis square
    end
end

val = 1;
M = VALS(:,:,val)';
y = repmat(1:15,10,1)'; % generate x-coordinates
x = repmat(1:10,15,1); % generate y-coordinates
% Generate Labels
t = num2cell(M); % extact values into cells
t = cellfun(@num2str, t, 'UniformOutput', false); % convert to string
% Draw Image and Label Pixels
figure
imagesc(M)
text(x(:), y(:), t, 'HorizontalAlignment', 'Center')