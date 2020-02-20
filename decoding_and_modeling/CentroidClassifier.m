%% CentroidClassifier
%
% Program to take in two sets of Centroids, undersample their neurons and
% then assign one set to their nearest neighbours in another set, we define
% neighbour as nearest angle now! We then output the accuracy and a measure
% of deviation.
%
% INPUT:
% y and z - two Neuron x Odour Centroid Representations
%   We use y to classify z!
% NumNeu - Number of Neurons to subsample
% Num - Number of iterations to run over for neuronal subsampling
%
% OUTPUT:
% Accuracy - average number of the odours correctly assigned
% Deviation - a measure of error
%
function [Accuracy, Deviation] = CentroidClassifier(y, z, NumNeu, Num)

% Setup
Neurons = size(y,1);
Odours = size(y,2);
Accuracies = zeros(1, Num);

% Iterate
for i = 1:Num
    % Subsample
    Subsample = randperm(Neurons, NumNeu);
    y_S = y(Subsample, :);
    z_S = z(Subsample, :);
    
    Assignment = zeros(1,15);
    % Assign each centroid
    for j = 1:Odours
        mindist = Inf;
        for k = 1:Odours
            dist = norm(z_S(:,j)-y_S(:,k));
            if dist < mindist
                mindist = dist;
                Assignment(j) = k;
            end
        end
        if Assignment(j) == j
            Assignment(j) = 1;
        else
            Assignment(j) = 0;
        end
    end
    
    Accuracies(i) = sum(Assignment)/Odours;
end

Accuracy = mean(Accuracies);
Deviation = std(Accuracies);