%% findCentroids.m
%
% Programme that takes in dataset and labels and outputs the centroids of
% each class
%
% INPUTS
%   Either
% X - dataset (N x P matrix)
% Y - labels (N x 1 vector)
%   Or
% X - dataset (C x T x P matrix)
%
% Where P is number of neurons, N is data points, C is classes, T is trials
% 
% OUTPUTS
%
% mu - centroids (C x P matrix)
%

function [mu] = findCentroids(X, Y)
    
if length(size(X)) == 2
    NumClasses = length(unique(Y));
    Classes = unique(Y);
    Trials = size(X, 1);
    Neurons = size(X, 2);
    mu = zeros(NumClasses, Neurons);

    for j = 1:NumClasses
        num = 0;
        for k = 1:Trials
            if Y(k) == Classes(j)
                mu(j, :) = mu(j,:) + X(k,:);
                num = num + 1;
            end
        end
        mu(j,:) = mu(j,:)/num;
    end
    
elseif length(size(X)) == 3
    NumClasses = size(X, 1);
    Repeats = size(X, 2);
    Neurons = size(X, 3);
    mu = zeros(NumClasses, Neurons);
    for j = 1:NumClasses
        for k = 1:Repeats
            mu(j, :) = mu(j,:) + squeeze(X(j,k,:))';
        end
        mu(j,:) = mu(j,:)/Repeats;
    end
end
end

