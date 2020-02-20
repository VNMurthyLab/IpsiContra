%% CentroidClassifierSide
%
% Program to take in two sets of Centroids, undersample their neurons and
% then run through the pairs of trials and assign them to the closest side
% centroid. Outputs accuracy
%
% INPUT:
% y and z - two Neuron x Odour Centroid Representations, order unimportant
% NumNeu - Number of Neurons to subsample
% Num - Number of iterations to run over for neuronal subsampling
%
% OUTPUT:
% Accuracy - average number of the odours correctly assigned
% Deviation - a measure of error
%
function [Accuracy, Deviation] = CentroidClassifySide(y, z, NumNeu, Num)

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
    
    Assignment = zeros(2, 15);
    % Assign each centroid
    for j = 1:Odours
        mindist = [Inf, Inf];
        dist = [0,0];
        classes = 1:15;
        classes(j) = [];
        
        CentroidIpsi = mean(y_S(:, classes),2);
        CentroidContra = mean(z_S(:,classes),2);
        Centroids = [CentroidIpsi, CentroidContra];
        for k = 1:2
            dist(1) = norm(y_S(:,j) - Centroids(:,k));
            dist(2) = norm(z_S(:,j) - Centroids(:,k));
            for m = 1:2
                if dist(m) < mindist(m)
                    mindist(m) = dist(m);
                    Assignment(m, j) = k;
                end
            end
        end
        for m = 1:2
            if Assignment(m, j) == m
                Assignment(m,j) = 1;
            else
                Assignment(m,j) = 0;
            end
        end 
    end
    
    Accuracies(i) = sum(Assignment(:))/(2*Odours);
end

Accuracy = mean(Accuracies);
Deviation = std(Accuracies);