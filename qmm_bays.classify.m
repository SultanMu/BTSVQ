function [posterior, class] = gmm_bays_classify (data, sD1,sM1, sD2,sM2)
% classifies given data into two classes based on the gaussian kernal parameters seleclted for both classes
% All the map units in both Maps are taken as the centres of gaussian kernals, 
% and the parameters of the mixture model are learned which are then used in the Basian frame work classification

%   data = matrix
          %the data to be classified in to two classes, should have the same dimentionality as of the two training
          %classes
%   sD1  = Data Struct 1
%   sM1  = Map Struct 1
%   sD2  = Data Struct 2
%   sM2  = Map Struct 2

%   Posteriors = Posterior probablities of the data for two classes
%   class = 0 and 1 for two classes


% Created by Mujahid Sultan Oct 19, 2001
% msultan@uhnres.utoronto.ca
%#######################################################

% find genes of training set ONE
% get the gene labels which are in the sM1
train1_in_data = sorting_general ('temp', sD1.labels,sM1.labels);
% now pick up the whole data vectors from the full data set on the above gene labels
train1_class = sD1.data(sM1.labels(train1_in_data),:);

% Repeat for the Training set TWO
train2_in_data = sorting_general ('temp',sD2.labels,sM2.labels);
train1_class = sD2.data(sM2.labels(train2_in_data);

% find the kernal parameters for both classes
[K1, P1]= som_estimate_gmm_2(train1_class,sM1,sD1);
[K2, P2]= som_estimate_gmm_2(train2_class,sM2,sD2);

% classify the data
[posterior, class]=gmm_em_prob_classifier_3b(K1,P1,K2,P2,train1_class,train2_class,data);
plot (posterior(:,1),',','FontSize','14')
