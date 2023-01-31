load('imports-85')
Y = X(:,14:15);
[n,d] = size(Y);
X1 = zscore(X(:,3));
X2 = zscore(X(:,7));
X3 = X(:,18)==20;
Xmat = [ones(n,1) X1 X2 X3];

Xcell = cell(1,n);
for i = 1:n
    Xcell{i} = [kron([Xmat(i,:)],eye(d))];
end

[beta,sigma,E,V] = mvregress(Xcell,Y);
beta