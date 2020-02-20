function code_for_figure_s3()

%%%%%%%%%%%%
%FIGURE S3D%
%%%%%%%%%%%%

%Generate artificial Poisson trials from actual data (Figure S3D)
load('tetrodeRecordings_OC_2s_notBaselineSubstracted.mat')
A_Poisson = {};
for mouse = 1:10;
    M = A{mouse};
    M = mean(M,4);
    AP = [];
    for t = 1:14
        AP(:,:,:,t) = poissrnd(M);
    end
    s = size(AP);
    for neuron = 1:s(1)
        for odor = 1:16
            for side = 1:s(3)
                for rep = 1:s(4)
                    AP(neuron,odor,side,rep) = AP(neuron,odor,side,rep) - mean(AP(neuron,16,side,:));
                end
            end
        end
    end
    A_Poisson{mouse} = AP;
end

%Calculate bilateral correlations in artificially-generated data (Figure S3D)
A = A_Poisson;
Data = {};
warning ('off','all')
for mouse = 1:10
    M = A{mouse};
    s = size(M);
    for side = 1:2
        [mouse side]
        PVAL = [];
        for neuron = 1:s(1);
            X = mean(squeeze(M(neuron,1:15,side,1:7)),2);
            Y = mean(squeeze(M(neuron,1:15,side,8:14)),2);
            lm = fitlm(X',Y','linear');
            lm = lm.Coefficients;
            lm = table2array(lm);
            pval = lm(2,4);
            PVAL = [PVAL ; pval];
        end
        Data{mouse,side} = PVAL;
    end
end
warning ('on','all')

%Calculate percentage bilat-corr neurons in artificially-generated data (Figure S3D)
alpha = 0.05;
Percent = zeros(10,2);
for side = 1:2
    for mouse = 1:10
        Pval = Data{mouse,side};
        isSignif = find(Pval<=alpha);
        Percent(mouse,side) = 100*length(isSignif)/length(Pval);
    end
end

%Plot percentage bilat-corr neurons in artificially-generated data (Figure S3D)
figure
mouse = [1 1 1 2 2 2 2 3 3 3];
colorPlot = {'r','b','k'};
for side = 1:2
    for region = 1:3
        subplot(1,2,side)
        hold on
        toPlot = Percent(find(mouse==region),side);
        plot(mouse(find(mouse==region)),toPlot,'o','MarkerEdgeColor',colorPlot{region})
        ylim([0 100])
    end
end
% Note: bootstrap analysis identical to Figure 4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%
%FIGURE S3B%
%%%%%%%%%%%%

Percentages = zeros(10,1);
alpha = 0.05;

%Percentage of bilaterally-correlated neurons
% Keep only neurons responding to at least 1 odor on 1 side (Figure S3B)
load('tetrodeRecordings_OC_2s.mat')
for mouse = 1:10;
    mouse
    PVAL = [];
    M = A{mouse};
    M = mean(M,4);
    for odor = 1:15
        M(:,odor,:) = M(:,odor,:) - M(:,16,:);
    end
    Mi = M(:,1:15,2);
    Mc = M(:,1:15,1);
    
    whichSignif = C{mouse};
    whichSignif = sum(sum(whichSignif,3),2);
    whichSignif = find(whichSignif>=1);
    Mi = Mi(whichSignif,:);
    Mc = Mc(whichSignif,:);
    
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
    [length(find(PVAL<=alpha))*100/length(PVAL) length(find(PVAL<=alpha)) length(PVAL)]
end

figure
plot([1 1 1 2 2 2 2 3 3 3],Percentages,'o')
ylim([0 50])
% Note: bootstrap analysis identical to Figure 4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%
%FIGURE S3C%
%%%%%%%%%%%%

Percentages = zeros(10,1);
alpha = 0.05;

%Percentage of bilaterally-correlated neurons
% Keep only neurons responding to at least 1 odor on each side (Figure S5Aiii)
load('tetrodeRecordings_OC_2s.mat')
for mouse = 1:10;
    mouse
    PVAL = [];
    M = A{mouse};
    M = mean(M,4);
    for odor = 1:15
        M(:,odor,:) = M(:,odor,:) - M(:,16,:);
    end
    Mi = M(:,1:15,2);
    Mc = M(:,1:15,1);
    
    whichSignif = C{mouse};
    whichSignifi = sum(whichSignif(:,:,2),2)&1;
    whichSignifc = sum(whichSignif(:,:,1),2)&1;
    whichSignif = find(whichSignifi & whichSignifc);
    Mi = Mi(whichSignif,:);
    Mc = Mc(whichSignif,:);
    
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
    [length(find(PVAL<=alpha))*100/length(PVAL) length(find(PVAL<=alpha)) length(PVAL)]
end

figure
plot([1 1 1 2 2 2 2 3 3 3],Percentages,'o')
ylim([0 50])
% Note: bootstrap analysis identical to Figure 4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%
%FIGURES S3F - S3H%
%%%%%%%%%%%%%%%%%%%

load('tetrodeRecordings_OC_waveform.mat')
step = [50:9:300];
nboot = 1000;
figure

%Spike width at half maximum amplitude, per region (Figures S3F and S3H)
mouse_region = [1 1 1 2 2 2 2 3 3 3];
for region = 1:3
    WVwidth = [];
    for mouse = find(mouse_region==region)
        WVwidth = [WVwidth;SpikeWidth{mouse}];
    end
    subplot(4,1,region)
    histogram(WVwidth,step,'DisplayStyle','stairs')
    [~,p] = HartigansDipSignifTest(WVwidth,nboot);
    title(['Hartigan s Dip Test, p = ',num2str(p)])
    xlim([0 300])
end

%Spike width at half max amplitude, bilaterally-corr neurons (Figure S3H)
neuron_id = [];
alpha = 0.05;
load('tetrodeRecordings_OC_2s.mat')
for mouse = 1:7;
    mouse
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
            neuron_id = [neuron_id;mouse neuron];
        end
    end
end

WVwidth = [];
for mouse = 1:7
    all = SpikeWidth{mouse};
    n = neuron_id(find(neuron_id(:,1)==mouse),2);
    WVwidth = [WVwidth;all(n)];
end
subplot(4,1,4)
histogram(WVwidth,step,'DisplayStyle','stairs')
[~,p] = HartigansDipSignifTest(WVwidth,nboot);
title(['Hartigan s Dip Test, p = ',num2str(p)])
xlim([0 300])

%Spike width at half max amplitude, box plots (Figure S3G)
cat = [];
val = [];
for mouse = 1:10
    M = BasalFR{mouse};
    val = [val;M];
    M = SpikeWidth{mouse};
    M(M<=100) = 0;
    M(M>100) = 1;
    cat = [cat;M];
end
for mouse = 1:7
    f = find(neuron_id(:,1)==mouse);
    f = neuron_id(f,2);
    M = BasalFR{mouse};
    M = M(f);
    val = [val;M];
    cat = [cat;2*ones(length(M),1)];
end
figure
boxplot(val,cat,'PlotStyle','compact')

[p,~,c] = kruskalwallis(val,cat);
multcompare(c)
