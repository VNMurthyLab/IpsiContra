function [R,P,A,B] = corrRPAB(X,Y)

s = size(X);
X = reshape(X,[s(1)*s(2) 1]);
Y = reshape(Y,[s(1)*s(2) 1]);

X(isnan(X)) = [];
Y(isnan(Y)) = [];

[r,p] = corrcoef(X,Y);
R = r(2,1);
P = p(2,1);

ab = polyfit(X,Y,1);
A = ab(1);
B = ab(2);

end
