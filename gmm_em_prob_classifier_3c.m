function [posteriors, posteriors_class] = gmm_em_prob_classifier_3c(K1,P1,K2,P2,sM1,sM2,X1,X2)

% Finds the posterior probablities of the codebook vectors X1 and X2 using the parameters of P1
% and P2. The parameters of k1 and k2 are generated by assuming one gaussian per
% data vector and finding its means and variences using EM algorithm. Then these parameters
% are used to classify a new vector.

% X1 = sM.codebook(:,[class1])
% X2 = sM.codebook(:,[class2])
% psteriors =  matrix of two columns 
%              the posteriors of the vector x interms of X1 and X2

% psteriors_class =  matrix of two columns 
%                  lists 1 for class X1 and 0 for class 2

[N1 NN1]= size(X1); 
[N2 NN2]= size(X2); 
R = (length(sM1.codebook) + length(sM2.codebook)); 
% Compute prior probability of each class 
PC1 = length(sM1.codebook)/R; 
%PC1 = .5; %as both classes are of same size
PC2 = 1-PC1; 
Ppr1 = P1;
Ppr2 = P2; 

[pD1,Pdm1,pmd1] = som_probability_gmm_2(X1, sM1, K1, P1);
[pD2,Pdm2,pmd2] = som_probability_gmm_2(X2, sM2, K2, P2);
% Compute posterior probabilities using Bayes Rule 

for ii = 1:N1 
	PC1X(ii) = pD1(ii)*PC1/(pD1(ii)*PC1+pD2(ii)*PC2); 
end 
for ii = 1:N2 
	PC2X(ii) = pD2(ii)*PC2/(pD1(ii)*PC1+pD2(ii)*PC2); 
end 

%PC2X = 1-PC1X; 
PCX = [PC1X',PC2X']; 
posteriors = PCX;

% Classify 
num = 0; 
for j = 1:length (PCX) 
	if PC1X(j) < 0.5 
		PC1X(j) = 0; 
		PC2X(j) = 1; 
		num = num + 1; 
	else 
		PC1X(j) = 1; 
		PC2X(j) = 0; 
	end 
end 
PCX_classs = [PC1X', PC2X'];
posteriors_class= PCX_classs;

% Print results
%Per1 = num/N*100 
% Plot results 
%for ii = 1:length(Mu) 
%makeEllipse(Mu(ii,1),Mu(ii,2),2*Sx(ii),2*Sy(ii)) 
%hold on 
%end 
%plot(X(:,1),X(:,2),'+ k') 
%xlabel('x','fontsize',16) 
%ylabel('y','fontsize',16) 
%print -depsc ClassTest.eps 