clear

load SimData.mat

X = data;

%%
figure(1); clf;
hist(X.rep,[1:10])