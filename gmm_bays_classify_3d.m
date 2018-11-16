function [labels1,labels2] = gmm_bays_classify_3d(sD0,sM0,data1,data2)
% function [labels1,labels2] = gmm_bays_classify_3d(sD0,sM0,data1,data2)

% This function generates two seperate mixture models taking the codebook 
% of the parent as centres of the voronoi regions and adding gausian noise to them.
% The PROCEDURE:
% 
% We take the codebook vector of the parent and generate a mixture model 1. the process
% of generating mixture model is to check the BMU for a given mi for all x. (where mi 
% are codebook vectors and x is a data vector under consideration.



% Data1= is the data though subset of sM0.codebook but only that indices which apear in
% sM11.labels (or sM11.codebook for that mater). Similarly the data 2 is for the other child.
%         
% the main idea here is that if data vector (or label) selectecd by the vector quantizer
% at a lowe level (at a child level) pridicted by the mixture model trained at the parent
% level (though using only the vectors in child, with full dimentionality), we can then deduce
% that the quantizer build on the smaller dimentionality data set and the partitive clustering
% (in this case partitive k-means) generated similar results. 
% 
% If the probabilities do not compromise to some optimal, iterate.
        
% sD, sM = data struct and trained map struct
% 
% Mujahid sultan, msultan@uhnres.utoronto.ca 
% beta 1.0   Feb 2003
%


% find the kernal parameters for both classes
[K1, P1]= som_estimate_gmm(sM0,data1);
[K2, P2]= som_estimate_gmm(sM0,data2);
%K2= K1; P2 = P1;

[N1 NN1]= size(data1); 
[N2 NN2]= size(data2); 

[pD1,Pdm1,pmd1] = som_probability_gmm(data1, sM0, K1, P1);
[pD2,Pdm2,pmd2] = som_probability_gmm(data2, sM0, K2, P2);
% Compute posterior probabilities using Bayes Rule 

prob1 = zeros(N1,1);
ind1 =zeros(N1,1);
for i = 1:N1
[prob1(i),ind1(i)] = max(pmd1(:,i));
end
labels1 = ind1(find(prob1>.50));
labels1 = sM0.labels(labels1,1)

prob2 = zeros(N2,1);
ind2 =zeros(N2,1);
for i = 1:N2
[prob2(i),ind2(i)] = max(pmd2(:,i));
end
labels2 =  ind2(find(prob2>.50));
labels2 = sM0.labels(labels2,1)

