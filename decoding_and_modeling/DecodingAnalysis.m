%% DecodingAnalysis.m
% This script generates all the decoding analysis for the paper:
% "Matched Bilateral Odor Responses in the Olfactory Cortex Point to 
% Non-Random Connectivity", Grimaud et al.
%
% First we import the data. 
% We then perform four decoding tasks:
% Task 1: Odour decoding in AON using different time periods, fig S5A
% Task 2: Odour decoding in cortical region and across hemisphers, fig 5B & C
% Task 3: Side decoding in all cortical regions, fig 5D
% Task 4: Odour decoding in each mouse and across hemipsheres, fig S5B

%% Import the Raw Data
% Tetrode recordings from cortex and OB
load Data_and_Notes/tetrode_recordings_formatted_for_alex_v2_ordered.mat
load Data_and_Notes/tetrode_recordings_formatted_for_alex_v2_ordered_ob.mat

% And all the times of sniffing
load Data_and_Notes/SniffTime.mat
load Data_and_Notes/SniffTimeOB.mat

% Collect them together for ease of passing around (just saved typing)
ArrayOB = [ArrayOB1; ArrayOB2; ArrayOB3];
DataSets = {ArrayOB, ArrayAON1, ArrayAON2, ArrayAON3, ArrayAPC1, ...
    ArrayAPC2, ArrayAPC3, ArrayAPC4, ArrayPPC1, ArrayPPC2, ArrayPPC3};

% Again just collecting SniffTimes together
OBSniffTime = cell(1, 375, 16, 3);
OBSniffTime(:, 1:27, :,:) = SniffTimeOB;
SniffTimeTot = cat(1,OBSniffTime, SniffTime);

% Neurons is a list of the number of neurons in each dataset element 
Neurons = zeros(1, length(DataSets));
for i = 1:length(DataSets)
    Neurons(i) = size(DataSets{i}, 1);
end

% A useful list to have lying around
Cortices = {'OB', 'AON', 'AON', 'AON', 'APC', 'APC', 'APC', 'APC',...
    'PPC', 'PPC', 'PPC'};

% And now just cleanup your workspace
clear ArrayAON1 ArrayAON2 ArrayAON3 ArrayAPC1 ArrayAPC2 ArrayAPC3 ...
    ArrayAPC4 ArrayPPC1 ArrayPPC2 ArrayPPC3 SniffTime SniffTimeOB OBSniffTime ...
    ArrayOB1 ArrayOB2 ArrayOB3

%% Now we start performing decoding - this case varying bin size (Fig S5A)
% We decode using the summed # of each neuron's spikes during a time bin 
% starting when the mouse sniffs. In this section we vary the size of the
% bin before seeing that 2 seconds is quite good.
%
% All decoding and data manipulation is done within the function
% NowClassify.m that accepts a list of parameter choices allowing you
% easily to switch the experiment being performed.
%
% In this instance we shall use it to find a good bin size.
%
% For meaning of these parameters see NowClassify.m (hopefully some are
% obvious though!)
Mice = 2:4;
Sides_Train = 3;
Sides_Test = 3;
Chosen_Odours = 1:15;
Num = 200;
Neuron = sum(Neurons(Mice));
Sniff_Align = 1;
Task = 0;

% The different Bin Sizes to run through
Bin_Sizes = 0.1:0.1:5;

% Create some empty arrays to store the results
Accuracies = zeros(1, length(Bin_Sizes)+1);
Deviations = Accuracies;
% This corresponds to 0 bin length, therefore chance (makes plot better)
Accuracies(1) = 1/length(Chosen_Odours);

% Now run through different bin sizes and test the odour classification acc
for k = 1:length(Bin_Sizes)
    Bin_Size = Bin_Sizes(k);
    disp(['Bin Size: ', num2str(Bin_Sizes(k))])
    [Accuracies(k+1), Deviations(k+1), ~, ~] = NowClassify(Mice, Sides_Train,...
        Sides_Test, Chosen_Odours, Bin_Size, Num, Neuron, Sniff_Align, Task, ...
        DataSets, SniffTimeTot);
end

% Plot the results
figure
errorbar([0,Bin_Sizes], Accuracies, Deviations)
ylabel('Accuracy')
xlabel('Bin Size')
title('Plot 1: Changing Spike Summing Bin effect on Accuracy (Fig S5A)')

% Finally save and plot up the results
save ('Data_and_Notes\ChangingBinSize', 'Accuracies', 'Deviations')

%% Next experiment: in each cortical area try decoding the odours. (Fig 5B & C)
% We combine all mice measurements from the same area into one
% psuedopopulation. We then train and decode across different neuronal
% populations, testing and training on ipsi and contra in all combinations.
%
% From here on it is very similar to the first bin size code.
Bin_Size = 2;
Cortices_Names = {'AON', 'APC', 'PPC'};
figure

% The Accuracies and deviations will be 3 cells, one corresponding to each
% cortical area. AON, APC, PPC
Accuracies = cell(3, 1);
Deviations = Accuracies;
% There will be twelve runs in total, 4 for each cortical area.
% 1 - contra train, contra test
% 2 - contra train, ispi test
% 3 - ipsi train, contra test
% 4 - ipsi train, ispi test
% For each we shall run across many neuronal populations
% These two will hold confusion matrices and the spread of accuracies at
% the full neuronal complement.
Confusions = zeros(3, 4, length(Chosen_Odours), length(Chosen_Odours)); 
AccSpread = zeros(3, 4, 200);

for Cortex = 1:3
    if Cortex == 1
        Mice = 2:4;
    elseif Cortex == 2
        Mice = 5:8;
    else
        Mice = 9:11;
    end
    Neuron = sum(Neurons(Mice));
    NeuronA = 10:10:50;
    NeuronB = 100:50:Neuron;
    Neuronal = [NeuronA, NeuronB, Neuron];      % Neuronal subsampling values

    ThisCortexAcc = zeros(4, length(Neuronal)+1);
    ThisCortexDev = ThisCortexAcc;
    ThisCortexAcc(:,1) = 1/(length(Chosen_Odours));

    for i = 1:4
        if i < 3
            Sides_Train = 2;
        else
            Sides_Train = 3;
        end
        Sides_Test = 3 - mod(i, 2);
    
        for k = 1:length(Neuronal)
            NumNeu = Neuronal(k);
            disp(['Cortex: ', Cortices{Mice(1)}, ', Side: ', num2str(i),', Neurons: ', num2str(Neuronal(k))])
            [ThisCortexAcc(i,k+1), ThisCortexDev(i,k+1), Confusion, Spread] = NowClassify(Mice, ...
                Sides_Train, Sides_Test, Chosen_Odours, Bin_Size, Num, NumNeu,...
                Sniff_Align, Task, DataSets, SniffTimeTot);
        
            % At one point we were comparing the spread of accuracies at a
            % particular number of neurons, I believe this is no longer in
            % the paper.
            if Neuronal(k) == 350
                AccSpread(Cortex, i,:) = Spread;
            end
            
            if k == length(Neuronal)
                Confusions(Cortex, i,:,:) = Confusion;
            end
            
        end
    end
    
    % ongoing plotting
    subplot(2, 3, Cortex)
    hold on
    errorbar([0, Neuronal], ThisCortexAcc(4,:), ThisCortexDev(4,:),'*', 'DisplayName','Test Ipsi')
    errorbar([0, Neuronal], ThisCortexAcc(3,:), ThisCortexDev(3,:),'*', 'DisplayName','Test Contra')
    title([Cortices_Names{Cortex}, ', Train: Ipsi']) 
    legend
    if Cortex == 1
        ylabel('Accuracy')
    end
    hold off
    
    subplot(2, 3, Cortex+3)
    hold on
    errorbar([0, Neuronal], ThisCortexAcc(1,:), ThisCortexDev(1,:),'*', 'DisplayName','Test Contra')
    errorbar([0, Neuronal], ThisCortexAcc(2,:), ThisCortexDev(2,:),'*', 'DisplayName','Test Ipsi')
    title([Cortices_Names{Cortex}, ', Train: Contra']) 
    legend
    if Cortex == 2
        xlabel('Number of Neurons')
    end
    if Cortex == 1
        ylabel('Accuracy')
    end
    hold off
    
    Accuracies{Cortex} = ThisCortexAcc;
    Deviations{Cortex} = ThisCortexDev;
end

sgtitle('Plot 2: Odour Decoding by Cortex, Fig 5B and C')
save ('Data_and_Notes\ByCortexOdourDecoding', 'Accuracies', 'Deviations', 'AccSpread', 'Confusions')

%% Third experiment: can we decode odour presentation side? (Fig 5D)
% Set task to 2, classifying the side
Task = 2;
figure
% And we'll only try and tell the difference between ipsi and contra
Sides_Train = [2,3];
Sides_Test = [2,3];

% Again the three cells are the three cortices
Accuracies = cell(3, 1);
Deviations = Accuracies;
Confusions = zeros(3, length(Sides_Train), length(Sides_Train)); 
AccSpread = zeros(3, 200);

for Cortex = 1:3
    if Cortex == 1
        Mice = 2:4;
    elseif Cortex == 2
        Mice = 5:8;
    else
        Mice = 9:11;
    end
    Neuron = sum(Neurons(Mice));
    NeuronA = 10:10:50;
    NeuronB = 100:50:Neuron;
    Neuronal = [NeuronA, NeuronB, Neuron];

    ThisCortexAccSide = zeros(1, length(Neuronal)+1);
    ThisCortexDevSide = ThisCortexAccSide;
    ThisCortexAccSide(1) = 1/(length(Sides_Train));

    for k = 1:length(Neuronal)
        NumNeu = Neuronal(k);
        disp(['Cortex: ', Cortices{Mice(1)}, ', Neurons: ', num2str(Neuronal(k))])
        [ThisCortexAccSide(k+1), ThisCortexDevSide(k+1), Confusion, Spread] = NowClassify(Mice, ...
            Sides_Train, Sides_Test, Chosen_Odours, Bin_Size, Num, NumNeu,...
            Sniff_Align, Task, DataSets, SniffTimeTot);

        if Neuronal(k) == 350
            AccSpread(Cortex,:) = Spread;
        end

        if k == length(Neuronal)
            Confusions(Cortex,:,:) = Confusion;
        end

    end
   
    % ongoing plotting
    subplot(1, 3, Cortex)
    hold on
    errorbar([0, Neuronal], ThisCortexAccSide, ThisCortexDevSide,'*')
    title(Cortices_Names{Cortex}) 
    if Cortex == 1
        ylabel('Accuracy')
    end
    if Cortex == 2
        xlabel('Number of Neurons')
    end
    hold off
    
    Accuracies{Cortex} = ThisCortexAccSide;
    Deviations{Cortex} = ThisCortexDevSide;
end

sgtitle('Plot 3: Side Decoding, Fig 5D')
    
save ('Data_and_Notes\SideDecoding', 'Accuracies', 'Deviations', 'AccSpread', 'Confusions')

%% Final Experiment: we look at the decoding in each individual mouse (Fig S5B)
% This is the same as Fig 5B and C except we decode using only one mouse, easy!
Task = 0;
figure

% There are now 10 storage cells, for each of the ten cortically measured mice
Accuracies = cell(10, 1);
Deviations = Accuracies;
% And again on full neuronal complement we store the four experiments' confusion
Confusions = zeros(10, 4, length(Chosen_Odours), length(Chosen_Odours)); 
AccSpread = zeros(10, 4, 200);

for Cortex = 2:11
    Mice = Cortex;
    Neuron = Neurons(Cortex);
    NeuronA = 10:10:50;
    NeuronB = 100:50:Neuron;
    Neuronal = [NeuronA, NeuronB, Neuron];

    ThisMouseAcc = zeros(4, length(Neuronal)+1);
    ThisMouseDev = ThisMouseAcc;
    ThisMouseAcc(:,1) = 1/(length(Chosen_Odours));

    for i = 1:4
        if i < 3
            Sides_Train = 2;
        else
            Sides_Train = 3;
        end
        Sides_Test = 3 - mod(i, 2);
    
        for k = 1:length(Neuronal)
            NumNeu = Neuronal(k);
            disp(['Mouse: ', num2str(Cortex), ', Cortex: ', Cortices{Cortex}, ', Side: ', num2str(i),', Neurons: ', num2str(Neuronal(k))])
            [ThisMouseAcc(i,k+1), ThisMouseDev(i,k+1), Confusion, Spread] = NowClassify(Mice, ...
                Sides_Train, Sides_Test, Chosen_Odours, Bin_Size, Num, NumNeu,...
                Sniff_Align, Task, DataSets, SniffTimeTot);
        
            % For the mice case we used to compare the spread at 50 neurons
            if Neuronal(k) == 50
                AccSpread(Cortex-1, i,:) = Spread;
            end
            
            if k == length(Neuronal)
                Confusions(Cortex-1, i,:,:) = Confusion;
            end
            
        end
    end
    
    % continual plotting
    subplot(5, 4, 2*Cortex - 3)
    hold on
    errorbar([0, Neuronal], ThisMouseAcc(4,:), ThisMouseDev(4,:),'*', 'DisplayName','Test Ipsi')
    errorbar([0, Neuronal], ThisMouseAcc(3,:), ThisMouseDev(3,:),'*', 'DisplayName','Test Contra')
    title([Cortices{Cortex}, ', Train: Ipsi']) 
    legend
    if mod(Cortex, 2) == 0
        ylabel('Accuracy')
    end
    if Cortex == 10 || Cortex == 11
        xlabel('Number of Neurons')
    end
    hold off
    
    subplot(5, 4, 2*Cortex - 2)
    hold on
    errorbar([0, Neuronal], ThisMouseAcc(1,:), ThisMouseDev(1,:),'*', 'DisplayName','Test Contra')
    errorbar([0, Neuronal], ThisMouseAcc(2,:), ThisMouseDev(2,:),'*', 'DisplayName','Test Ipsi')
    title([Cortices{Cortex}, ', Train: Contra']) 
    legend
    if Cortex == 10 || Cortex == 11
        xlabel('Number of Neurons')
    end
    hold off
    
    Accuracies{Cortex-1} = ThisMouseAcc;
    Deviations{Cortex-1} = ThisMouseDev;
end

sgtitle('Plot 4: Odour Decoding by Mouse, Fig S5B')

save ('Data_and_Notes\ByMouseOdourDecoding', 'Accuracies', 'Deviations', 'AccSpread', 'Confusions')