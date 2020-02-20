%% NowClassify.m
% A function for preparing our feature vectors then attempting to decode
% from them. Relatively easy to customise to perform different tests.
%
% One improvement that could be made is saving the same dataset that we
% keep neuronally subsampling from, never quite got around to it though!
%
% Parameter Input choices are as follows:
%
% INPUTS
%
% Mice: Which mice to be used in this experiment. (OB treated as one)
%    1 OB, 2 - 4 AON, 5 - 8 APC, 9 - 11 PPC
%
% Side_Train and Side_Test: Which nostrils odours are presented for
% training and testing the classifier respectively.
%   1 - both nostrils presented, 2 - contra, 3 - ipsi
%
% Chosen_Odours: set of odours to use.
%   Numbers from 1 to 15, as listed in the paper fig S1A.
%
% Bin_Size: Size of the bin to use.
%   Measured in seconds and 0 is the sniff time.
%
% Num: the number of iterations to run through.
%   To average over different train/test set or neuron subset choices.
%
% NumNeu: the number of neurons to use.
%
% SniffAlign: Whether to align 0s to the sniff or not.
%   1 - align, 0 - 0s = odour onset (not used in paper).
%
% Task: controls which task you would like to be attempted.
%   0: Just classify to a particular odour
%   1: Classify odour AND side
%   2: Classify just the side.
%
% Also takes as input all the tetrode data and all the sniff times.
%
% OUTPUTS
%
% The mean and deviation of the accuracy over all the Num iterations.
% The confusion matrix of the classification task. Only calculated when you
% are using the full complement of neurons (to save repeated slow calc).
% The full list of Num Accuracies (believe it is no longer used in paper)

function [Accuracy, Deviation, Confusion, Accuracies] = NowClassify(Mice,...
    Sides_Train, Sides_Test, Chosen_Odours, Bin_Size, Num, NumNeu, ...
    Sniff_Align, Task, DataSets, SniffTime)

% Will hold training and test data (X) and labels (Y)
X_Train = [];
Y_Train = [];
X_Test = [];
Y_Test = [];

% If you are trying to decode odour presentation side then you must train 
% on every side that you test on. Here we test for that.
if Task
    if length(Sides_Train) ~= length(Sides_Test)
        disp('Cannot classify sides without training on them')
        return
    end
    for i = 1:length(Sides_Train)
        counter = 0;
        for j = 1:length(Sides_Test)
            if Sides_Train(i) ~= Sides_Test(j)
                counter = counter + 1;
            end
        end
        if counter == length(Sides_Test)
            disp('Cannot classify sides without training on them')
            return
        end
    end
end

% These two will ensure that in multi-side experiments (i.e. Task = 1 or 2)
% different sides will get different labels
counter_Train = 0;
counter_Test = 0;
% This will mark which test datasets have already been loaded while
% creating the training dataset, in order to avoid wasted effort.
Tester_Marker = [];

% Now extract the data required
for i = 1:length(Sides_Train)
    [X, Y] = DataCreator(Mice, Sides_Train(i), Chosen_Odours, Bin_Size,...
        Sniff_Align, DataSets, SniffTime);
    
    % Creating correct labels depending on the task
    if Task == 1
        Y2 = Y + 15*counter_Train;
        counter_Train = counter_Train+1;
    elseif Task == 2
        counter_Train = counter_Train+1;
        Y2 = repelem(counter_Train, 15);
    else
        Y2 = Y;
    end
    
    X_Train = [X_Train, X];
    Y_Train = [Y_Train, Y2];
    
    % Check if this dataset is also being used for testing, if so add.
    for j = 1:length(Sides_Test)
        if Sides_Train(i) == Sides_Test(j)
            Tester_Marker = [Tester_Marker, Sides_Test(j)];
            
            if Task == 1
                Y2 = Y + 15*counter_Test;
                counter_Test = counter_Test+1;
            elseif Task == 2
                counter_Test = counter_Test + 1;
                Y2 = repelem(counter_Test, 15);
            else
                Y2 = Y;
            end
            
            X_Test = [X_Test, X];
            Y_Test = [Y_Test, Y2];
        end
    end
end

% Now collect all the remaining test datasets (same process)
for i = 1:length(Sides_Test)
    flag = 0;   % marks presence in Tester_Marker
    if ~isempty(Tester_Marker)
        for k = 1:length(Tester_Marker)
            if Tester_Marker(k) == Sides_Test(i)
                flag = 1;
            end
        end
    end
    if ~flag
        [X, Y] = DataCreator(Mice, Sides_Test(i), Chosen_Odours, Bin_Size,...
            Sniff_Align, DataSets, SniffTime);

        if Task == 1
            Y = Y + 15*counter_Test;
            counter_Test = counter_Test+1;
        elseif Task == 2
            counter_Test = counter_Test + 1;
            Y = repelem(counter_Test, 15);
        end

        X_Test = [X_Test, X];
        Y_Test = [Y_Test, Y];
    end
end

% Work out the labels in this particular task
if Task == 0
    Labels = Chosen_Odours;
elseif Task == 1
    if length(Chosen_Odours) ~= 15
        disp('Not implemented: classifying both odours and sides for less than 15 odours')
    end
    Labels = 1:15*length(Sides_Train);
elseif Task == 2
    Labels = 1:length(Sides_Train);
end
NumTest = length(Labels);
Confusion = zeros(NumTest);     % Will hold confusion matrix

% If we are using all the neurons in a set there is no need to do multiple
% iterations over different neuronal populations, therefore do only Num_Neu_Itr.
Neurons = size(X_Train, 1);     % Total number of neurons, from which sampling
if NumNeu == Neurons
    Num_Neu_Itr = 1;
else
    Num_Neu_Itr = Num;
end

Accuracies = zeros(1, Num_Neu_Itr);    % will hold the goodies
Deviations = Accuracies;

% Run through each iteration choosing a different population of neurons,
% then pass to linear fitting function
for i = 1:Num_Neu_Itr
    Chosen_Neurons = randperm (Neurons, NumNeu);
    X_Train_NSub = X_Train(Chosen_Neurons,:,:);
    X_Test_NSub = X_Test(Chosen_Neurons,:,:);

    if Task < 2
        [Accuracies(i), Deviations(i), assignments] = Centroid_Leave_Out_One_Per_Odour(X_Train_NSub, Y_Train, X_Test_NSub, Y_Test);
    else
        [Accuracies(i), Deviations(i), assignments] = Centroid_Leave_Out_One_Per_Odour_Side(X_Train_NSub, Y_Train, X_Test_NSub, Y_Test);
    end
    
    % Use the assignments to calculate the confusion matrix
    if NumNeu == Neurons
        for n = 1:size(assignments, 1)
            index1 = find(Chosen_Odours == assignments(n,1));
            index2 = find(Chosen_Odours == assignments(n,2));
            Confusion(index1,index2) = Confusion(index1,index2) + 1;
        end
        Normaliser = sum(Confusion(1,:));
        Confusion = Confusion/Normaliser;
    end
end

% Finally calculate our average accuracy
Accuracy = mean(Accuracies);
% Combine all the deviations across different train/test sets in quadrature
% sometimes this is meaningless ( = 0 since no variation from train and
% test set choice) then we just use the deviation over the accuracies.
Deviation = (mean(Deviations.^2))^0.5;
if Deviation == 0
    Deviation = std(Accuracies);
end

disp(['Accuracy = ', num2str(Accuracy),' +/- ', num2str(Deviation), ' (Chance = ', num2str(1/NumTest), ')'])
end
