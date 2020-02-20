%% Centroid_Leave_Out_One_Per_Odour.m
%
% Here we take in a set of training data with labels, regroup it into
% trials, we then remove one of the data points from each class for training, 
% before the same removed points but from the test data sets.
%
% The classifier at work groups the classes into centroids and assigns to
% the closest euclidean centroid.
%
% INPUTS
%
% X Datasets - Neurons x Odours x Trials
% Y Labels - Odours

function[Accuracy, Deviation, Assignments] = Centroid_Leave_Out_One_Per_Odour(X_Train, Y_Train, X_Test, Y_Test)

% Lets just get a collection of useful stuff ready
Classes = unique(Y_Train);
NumClasses = length(Classes);
if length(unique(Y_Test)) ~= NumClasses
    disp('Different numbers of odours in test and train! No can do!')
    return
end
Trials = size(X_Train, 3);
Repeats = size(X_Train, 2)/NumClasses;
Neurons = size(X_Train, 1);

% Now we sort out the data into the nice 2D form
X_Train = permute(X_Train, [3, 2, 1]);
X_Train = reshape(X_Train, [NumClasses*Repeats*Trials, Neurons]);
Train_Trials = size(X_Train, 1);
Y_Train = repelem(Y_Train, Trials);
X_Test = permute(X_Test, [3, 2, 1]);
X_Test = reshape(X_Test, [NumClasses*Repeats*Trials, Neurons]);
Y_Test = repelem(Y_Test, Trials);

Assignments = zeros(200*NumClasses, 2);
Accuracies = zeros(1, 200);
y_test = zeros(1, NumClasses);
x_test = zeros(NumClasses, size(X_Test, 2));

% Iterate through removing one different trial each time
for  i = 1:200
    L = 0;
    index1 = 0;
    index2 = 0;
    x_test = [];
    x_train = X_Train;
    y_train = Y_Train;
    y_test = [];
    removers = zeros(1, NumClasses);
    
    % Remove one bit of data from each odour class and test on it
    for k = 1:NumClasses
        remove = randperm(Trials*Repeats, 1);
        sort = find(Y_Test == Classes(k));
        remove2 = sort(remove);
        removers(k) = remove2;
        x_test(k, :) = X_Test(remove2,:);
        y_test(k) = Y_Test(remove2);
    end
    
    x_train(removers,:) = [];
    y_train(removers) = [];
    
    % Find the means of each class
    mu = zeros(NumClasses, Neurons);
    for j = 1:NumClasses
        num = 0;
        for k = 1:(Train_Trials-length(y_test))
            if y_train(k) == Classes(j)
                mu(j, :) = mu(j,:) + x_train(k,:);
                num = num + 1;
            end
        end
        mu(j,:) = mu(j,:)/num;
    end
    
    for k = 1:length(y_test)
        % Now assign the left over points to closest mean
        winner = Inf;
        assignment = 0;
        for j = 1:NumClasses
            dist = norm(x_test(k,:) - mu(j,:));
            if dist < winner
                winner = dist;
                assignment = Classes(j);
            end
        end

        % Before finally seeing if the assigment is correct
        if assignment == y_test(k)
            L = L + 1;
        end
        Assignments(k + (i-1)*NumClasses, :) = [y_test(k), assignment];
    end
    L = L/length(y_test);
    Accuracies(i) = L;
end

Deviation = std(Accuracies);
Accuracy = mean(Accuracies);
end