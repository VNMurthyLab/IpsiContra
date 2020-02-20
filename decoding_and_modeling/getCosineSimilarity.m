%% getCosineSimilarity.m
% program to calculate the cosine similarity of members of a matrix

function [CosineSims] = getCosineSimilarity(x)
N = size(x, 2);

% Normalise each response vector
for i = 1:N
    x(:,i) = x(:,i)/norm(x(:,i));
end

% Calculate the dot products
dotprod = x'*x;

% retrun N*(N-1)/2 vector
mask = tril(true(N),-1);
CosineSims = dotprod(mask);