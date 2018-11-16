function [level] = p_tree_classify_supervised (sD,labels, filename,print)
% Same as p_tree function, the addition is that it makes SOM at each level and saves sM0 at each level and sM and sD 
% fro every child of a node.
% The output is passed to ??????????????????????

% level_classify  = p_tree_classify(sDt, 'labels', 'filename.txt',0 or 1)
% level = data structure containing k-means results
% ith level contains 2^i children
% labels field is empty string (for future use)
% Level(1).child(1).data 
% Level(2).child(1).data
% 
% sD, sM = data struct and trained map struct
% dim = number of columens in the data set
% dlenght = number of rows in the data set
%

% Mujahid sultan, msultan@uhnres.utoronto.ca 
% beta 1.0
%############################################################################


% if isstruct(sD) 
%  data = sD.name; 
%  labels = sD.labels;
%  comp_names = D.comp_names;
%else 
%  data = inputname(1);
%  labels = inputname(2);
%end
%super_error = 0 , %initialize the supervised error
error1_old = 0;
error2_old = 0;
counter = 0;
while 1 %super_error  % the supervised condition
    
    % Divide the data set into two children using k-means
    % [centroids,clusters,err] = kmeans(method, D, n, [epochs], [verbose])
    % centroids are in the 'dim' dimention
    
    % make the SOM at root level
    [centroids, clusters,err] = kmeans ('batch',sD,2,'epochs',10000);
    clusters = [(1:length(sD.data(:,1)))' clusters]; % annotate clusters with indix
    
    
    % Find out the indices of the vectors (rows) for both classes 
    % Note that you can do it quite easily but this is the initial implementation, so I keep it
    l= 0;
    k= 0;

    clear A B;
    M = clusters(:,1);   
    for ii = 1:length (M)
        % check to which cluster it belongs, 1st or 2nd 
        if clusters(ii,2) == 1;
            l= l+1;
            % record the gene number for 1st class
            A([l]) = clusters(ii,1);
        elseif clusters(ii,2) == 2
            k = k +1;
            % record the gene number for 2nd class
            B([k]) = clusters(ii,1);
        end  
    end
    % make a data struct to store two children   
    
    
    clear level;
    sD0  =som_data_struct (sD.data','comp_names',sD.labels,'labels',sD.comp_names);
    sD1 = som_data_struct (sD.data(A,:)','comp_names',sD.labels(A),'labels',sD.comp_names);
    sD2 = som_data_struct (sD.data(B,:)','comp_names',sD.labels(B),'labels',sD.comp_names);
    
    sM0 = som_make(sD0);
    sM1 = som_make(sD1);%,'msize',sM.topol);
    sM2 = som_make(sD2);%,'msize',sM.topol);
    
    %@@@@@@@@@@@@@@@@@ these are the additions made to p_tree_classify
% % %     sM00 = som_autolabel_rank (sM0,sD0,'add1');
% % %     sM11 = som_autolabel_rank (sM1,sD1,'add1');
% % %     sM12 = som_autolabel_rank (sM2,sD2,'add1');         
% % %     
    %%%[inds11,inds12]= label_comp(sM00.labels(:,1),sM11.labels(:,1),sM12.labels(:,1)); 
    %############## first scenario
    % the question here is whether we should pass the cated_ranked_genes11 and 12 or the whole codebook for that partition of the data
    % in the origional implementaion we are passing the whole subset of the codebook of paretn, but see what happens when we pass only
    % cated)ranked ones
    
    [posterior, class1_labels,class2_labels] = gmm_bays_classify_3c(sM0.labels(:,1),sM0.codebook(:,A),sM0.codebook(:,B),sD1,sM1, sD2,sM2);   
    
    error1_new = class1_labels
    error2_new = class2_labels
    
    if error1_new == 0 | error2_new == 0, counter = counter+1,  end;
    if counter > 5; break; end;
    
    if error1_new >= error1_old & abs(error1_new - error1_old) <= .05
        break
    elseif error2_new >= error2_old & abs(error2_new - error2_old) <= .05
        break
    else
        error1_old = class1_labels;
        error2_old  = class2_labels;
    end
                        
end  % while
% % %          % now- the posterior's first column gives the probablity of this gene to belong to class1
% % %          % and 2nd column gives the prob of this gene to belong to the class2.
% % %          % how ca we now check the correctness- the answer is to compare the labels.
% % %          % but we need to go into the functions of probablities and check the indices of these probablities....to which these are associated
% % %          correct_label1 = sM00.labels(inds11);
% % %          predicted_label1 = class1_labels;
% % %          classification1  = strcmp(correct_label1,predicted_label1);
% % %          correct_classification1  = lenght(classification1(classification1==1))/length(classification1);
% % %          
% % %          correct_label2 = sM00.labels(inds12);
% % %          predicted_label2 = class2_labels;
% % %          classification2  = strcmp(correct_label2,predicted_label2);
% % %          correct_classification2  = lenght(classification2==1)/length(classification2);
% % %          %################## end first Secnario
% % %          
% % %          %################## Second Scnario
% % %          [p_cb1,k_cb1]  =som_estimate_gmm (sM00,data1);
% % %          [p_cb2,k_cb2]  =som_estimate_gmm (sM00,data2);
% % %          
% % %          
% % %          
% % %          
         
         
     %level(0).sD0 = sD0;
     %level(0).sM0 = sM0;
     level(1).child(1).data= sD.data(A,:);
     level(1).child(1).labels= [A];%sD.labels(A); % store vector(rows) labels as well
     level(1).child(1).label_names= sD.labels(A);
     level(1).child(1).sD = sD1;
     level(1).child(1).sM = sM1;

     level(1).child(2).data = sD.data(B,:);
     level(1).child(2).labels= [B];%sD.labels(A); % store vector(rows) labels as well
     level(1).child(2).label_names= sD.labels(B);     
     level(1).child(2).sD = sD2;
     level(1).child(2).sM = sM2;

     
     % print the labels of these two children on file
     if print == 1
     fid = fopen (filename,'w');
     fprintf (fid,'\n%s\n   >> ','level(1).child(1)--->level(2).child(1 & 2)');%int2str((counter*10)+1));
     fprintf(fid, '\n');
     for iii = A, fprintf (fid, '%s,  ',sD.labels{iii});end 
     fprintf(fid, '\n');
     %fclose (fid);
     %fid = fopen (filename,'w');
     fprintf (fid,'\n%s\n   >> ','level(1).child(2)--->level(2).child(2 & 3)');%int2str((counter*10)+1));
     fprintf(fid, '\n'); 
     for iii = B, fprintf (fid, '%s,  ',sD.labels{iii});end 
     fprintf(fid, '\n');
     fclose (fid);
     end

%level(0).child(1)= struct('data',sD.data);

i =1; % initialize the level loop
while 1 % loop till the breaking condition   
%for i = 1 : 8 % 50 levels of Partative tree (may be take it as argument at later stage)     
    try, % exception for the level error
        if isstruct(level(i))
            %fprintf ( 'Level \n',i)
            clear children;
            cc = 1; % Initialize the child return
            for j = 1:2^i
        
                %select the data for this partition
                %level(i)= struct('data',level(i).child(m).data);
            
                %fprintf(2,' level / child [%d, %d] \n',i, j);   
                %level_ = i
                %child_ = j
                %pause % press a key to continue
            try, 
                %if isstruct(level(i).child(j))
                if (level(i).child(j).data ~=[])%&(isstruct(level(i).child(j)))    
                    %s = length(level(i).child(j).data);
                
                    fprintf(1,' level / child [%d, %d] \n',i, j);   
                    children([cc]) = j;
                    cc = cc+1;
                    
                    error1_old = 0;
                    error2_old  = 0;
                    counter = 0;
                    while 1
                        
                        [centroids, clusters,err] = kmeans ('batch',level(i).child(j).data,2, 'epochs',10000);
                        clusters = [(1:length(level(i).child(j).data(:,1)))' clusters]; % annotate clusters with indix
                        % Find out the indices of the vectors (rows) for both classes 
                        l= 0;
                        k= 0;
                        
                        clear A B;
                        M = clusters(:,1);   
                        for ii = 1:length (M)
                            % check to which cluster it belongs, 1st or 2nd 
                            if clusters(ii,2) == 1;
                                l= l+1;
                                % record the gene number for 1st class
                                A([l]) = clusters(ii,1);
                            elseif clusters(ii,2) == 2
                                k = k +1;
                                % record the gene number for 2nd class
                                B([k]) = clusters(ii,1);
                            end  
                        end
                        
                        %****************** BREAKING CONDITION *********************************
                        % calculate STD of parent
                        std_parent = mean(std(level(i).child(j).data)); % take mean of parent (try different criterian, as well)
                        
                        % 1) 
                        %varDiff = abs(var (centroids(1,:)) - var(centroids(2,:))); % not good results
                        % comments ---- feels like it is not true representative of the children, stops too early
                        
                        % 2)
                        %varDiff = mean(var (level(i).child(j).data(A,:))-var(level(i).child(j).data(B,:)));
                        % comment ---- 
                        %----varDiff = abs(mean(var (level(i).child(j).data(A,:))-var(level(i).child(j).data(B,:))));
                        % comment ----- I think the best (on lung8)
                        
                        % 3)
                        %varDiff = mean(min(var (level(i).child(j).data(A,:)),var(level(i).child(j).data(B,:))));
                        %varDiff = mean(max(var (level(i).child(j).data(A,:)),var(level(i).child(j).data(B,:))));
                        % comment ---- MAX could not produce any result
                        %comments --- gave even bad results than 2
                        
                        % 4)
                        % This is the expresstion used in pat_cluster_intersection
                        % varDiff=  ( mean(varience(1,:))+  mean(varience(2,:)))/2
                        varDiff = (mean(var (level(i).child(j).data(A,:)))+mean(var(level(i).child(j).data(B,:))))/2;
                        %***************** END BREAKING CONDITION ********************************
                        
                        % break if the STD of parent becomes less that Difference of varience of the children.
                        %if (std_parent > varDiff) 
                            
                            
                            sD0 = som_data_struct (level(i).child(j).data','comp_names',level(i).child(j).label_names,'labels',sD.comp_names);                 
                            sD1 = som_data_struct (level(i).child(j).data(A,:)','comp_names',level(i).child(j).label_names(A),'labels',sD.comp_names);
                            sD2 = som_data_struct (level(i).child(j).data(B,:)','comp_names',level(i).child(j).label_names(B),'labels',sD.comp_names);
                            sM0 = som_make (sD0); 
                            sM1 = som_make(sD1);%,'msize',sM.topol);
                            sM2 = som_make(sD2);%,'msize',sM.topol);
                            
                            
                            [posterior, class1_labels,class2_labels] = gmm_bays_classify_3c(sM0.labels(:,1),sM0.codebook(:,A),sM0.codebook(:,B),sD1,sM1, sD2,sM2);   
                            
                            error1_new = class1_labels
                            error2_new = class2_labels
                            
                            
                            if error1_new == 0 | error2_new == 0, counter = counter+1,  end;
                            if counter > 5; break; end;
                            
                            if error1_new >= error1_old & abs(error1_new - error1_old) <= .05
                                break
                            elseif error2_new >= error2_old & abs(error2_new - error2_old) <= .05
                                break
                            else
                                error1_old = class1_labels;
                                error2_old  = class2_labels;
                            end
                            
                        end % while 
                       
                       level(i).child(j).sD0 = sD0;
                       level(i).child(j).sM0 = sM0;
                       level(i).child(j).labelsA = A;
                       level(i).child(j).labelsB = B;
                       level(i+1).child(2*j-1).label_names = level(i).child(j).label_names(A);
                       level(i+1).child(2*j-1).data = level(i).child(j).data(A,:);
           	           level(i+1).child(2*j-1).labels= level(i).child(j).labels(A);  
           	           level(i+1).child(2*j-1).sD = sD1;
                       level(i+1).child(2*j-1).sM = sM1;  
                       
                       level(i+1).child(2*j).label_names=  level(i).child(j).label_names(B);
                       level(i+1).child(2*j).data = level(i).child(j).data(B,:);
                       level(i+1).child(2*j).labels= level(i).child(j).labels(B);  
                       level(i+1).child(2*j).sD= sD2;  
                       level(i+1).child(2*j).sM= sM2;                       
                       
                     if print == 1 
                     fid = fopen (filename,'a');
                     fprintf (fid,'\n%s\n   >> ',(['level(' int2str(i+1) ').child(' int2str(2*j-1) ') ---> level(' ...
                                                            int2str(i+2) ').child(' int2str(4*j-3) ' & ' int2str(4*j-2) ')']));
                     fprintf(fid, '\n');                               
                     for iii = A, fprintf (fid, '%s,   ',level(i).child(j).label_names{iii});end 
                     fprintf(fid, '\n');
                     fclose (fid);

                     fid = fopen (filename,'a');
                     fprintf (fid,'\n%s\n   >> ',(['level(' int2str(i+1) ').child(' int2str(2*j) ') ---> level(' ...
                                                            int2str(i+2) ').child(' int2str(4*j-1) ' & ' int2str(4*j) ')']));
                     fprintf(fid, '\n');
                     for iii = B, fprintf (fid, '%s,   ',level(i).child(j).label_names{iii});end 
                     fprintf(fid, '\n');
                     fclose (fid);
                     end
                 
% %                     else
% %                      
% %                      level(i+1).child(2*j-1).data = [];
% %                      level(i+1).child(2*j).data =   [];
% %                         Level(i+1).child(2*j-1).lables = [];
% %                         Level(i+1).child(2*j).labels = [];    
% %                     end
            
                else
                    % go and try next j, if all j th children are empty then exit
                    fff= 0; 
                    %level(i+1).child(2*j-1).data = [];
   	                %level(i+1).child(2*j).data =   [];
                end %if 
 
            catch, % if the i exceeds the generated level
                %disp (lasterr) Never display this error as this error is wanted to go to next step
                fff = 0; % simply take next j value
                %break % break the for loop
            end % try
        end % for j
        level(i).children = children;
        i =i+1; 
%         if i >100
%             return
%         end
    else 
        break
    end % if
   
catch, % if the i exceeds the generated level
    return % break the for loop and return the calling function
end % try
end % for or while      