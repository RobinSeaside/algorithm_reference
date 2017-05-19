close all; clear;
d = 2;
k = 3;
n = 500;
[X,label] = kmeansRnd(d,k,n);
y = kmeans(X,k);
plotClass(X,label);
figure;
plotClass(X,y);

%%
close all; clear;
tmp = load('K:/NTU_thesis/MyAlgorithms/code_reference/kmeans/kmeans/cache/x.mat');
X = tmp.test;
k = 2;
y = kmeans(X,k);
figure;
plotClass(X,y);