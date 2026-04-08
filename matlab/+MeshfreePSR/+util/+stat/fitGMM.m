function [GMM]=fitGMM(X,maxGauss,reg_GMM)
GMModels = cell(1,maxGauss);
options = statset('MaxIter',1500,'Display','off');
AIC=zeros(maxGauss,1); %this is the AIC criterion which is used to select the number of Gaussians
seed=73; 
rng(seed);%for reproducibility (wassGMM is deterministic but fitgmdist is not...)
for k=1:maxGauss
    GMModels{k} = fitgmdist(X,k,'Options',options,'RegularizationValue',reg_GMM);
    AIC(k)= GMModels{k}.AIC;
end
[~,numComponents] = min(AIC);
GMM=GMModels{numComponents};  
end