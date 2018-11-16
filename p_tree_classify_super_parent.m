function [level,level0] = p_tree_classify_super_parent(sD,labels, filename,print)
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
% beta 1.0   Feb 2003
%


% if isstruct(sD) 
%  data = sD.name; 
%  labels = sD.labels;
%  comp_names = D.comp_names;
%else 
%  data = inputname(1);
%  labels = inputname(2);
%end

%^^^^^^^^^^ROOT LEVEL^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%super_error = 0 , %initialize the supervised error
prod_old = 100*100/(100+100);
prod_max = 0;
counter = 0; counter0 = 0; counterNan = 0; max_counter = 0; iter = 1;
%fprintf(1,' level / child [%d, %d] \n',il+1, jl);   
while 1 % the supervised condition, iterate till the conditions are satisfied
    % in every iteratioon partition the data with k-means
    
    % Divide the data set into two children using k-means
    % centroids are in the 'dim' dimention
    [centroids, clusters,err] = kmeans ('batch',sD,2,'epochs',10000);
    clusters = [(1:length(sD.data(:,1)))' clusters]; % annotate clusters with indix
    
    % Find out the indices of the vectors (rows) for both classes 
    % Note that you can do it quite easily but this is the initial implementation, so I keep it
    l= 0; k= 0; clear A; clear B; 
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
    
    % Changes made in Nov 2003
    
    sM0 = som_make(sD0);
    sM0 = checkTopol (sM0);
    toPology = checkTopol (sM0);
    sM1 = som_make(sD1);%,'msize',sM.topol);
    sM1 = checkTopol (sM1);
    sM2 = som_make(sD2);%,'msize',sM.topol);
    sM2 = checkTopol (sM2);

    %@@@@@@@@@@@@@@@@@ these are the additions made to p_tree_classify in origional p-tree
    % Changes made in Nov 2003, to make run fast...this is origional                                
%                                 sM00 = som_autolabel_rank (sM0,sD0,'add1');
%                                 sM11 = som_autolabel_rank (sM1,sD1,'add1');
%                                 sM12 = som_autolabel_rank (sM2,sD2,'add1');         

    sM00 = som_autolabel (sM0,sD0,'auto');
    sM11 = som_autolabel (sM1,sD1,'auto');
    sM12 = som_autolabel (sM2,sD2,'auto');
    
    
    %     % SCNARIO 1) find the intersections and their indices in the parent codebook sM00 is the parent
    %     [inds11,inds12]= label_comp(sM00.labels(:,1),sM11.labels(:,1),sM12.labels(:,1)); 
    %     % find the probabilities of pridiction
    %     [class_labels1,class_labels2] = gmm_bays_classify_3d(sD0,sM0,sM0.codebook(inds11,:),sM0.codebook(inds12,:));   
    
    % SCNARIO 2) find the intersections and their indices in the parent all
    % data 
    [inds11,inds12]= label_comp(sD0.labels,sM11.labels(:,1),sM12.labels(:,1)); 
    % find the probabilities of pridiction
    [class_labels1,class_labels2] = gmm_bays_classify_3d(sD0,sM00,sD0.data(inds11,:),sD0.data(inds12,:));   
    
    
    
    % the scores are the actual labels and the pridicted labels
    score1 = 100*length(intersect(class_labels1,sD0.labels(inds11)))/length(inds11)
    score2 = 100*length(intersect(class_labels2,sD0.labels(inds12)))/length(inds12)
    
    if score1 >= 95 & score2 >= 95; level0.child0.prod(iter) = 50; break;  end;
    if isnan(score1)==1 | isnan(score2)==1, counterNan= counterNan+1
        if counterNan > 3
            if exist('optimal_A')==1 && exist('optimal_B')==1
            A = optimal_A;
            B = optimal_B; 
           end
            break; 
        end 
    end
    if score1==0 | score2==0, counter0= counter0+1
        if counter0 > 3 && exist('optimal_A')==1 && exist('optimal_B')==1
            A = optimal_A;
            B = optimal_B; 
            break; 
        end
    end
    if score1==0, score1=1; elseif score2==0, score2=1;end;
    
    % prod is the product of both scores divided by the sum
    % Maximizing prod will garuntee optmal solution
    if isnan(score1)==0 & isnan(score2)==0
        
        prod_new = (score1 * score2)/(score1 + score2) % maximum value of prod can take is 50
        
        if prod_new >= 40,level0.child0.prod(iter) = prod_new; break; end;                                                             
        
        % set the max_prod to maximum of prod_old and prod_new
        if prod_new <= 40 
            if prod_new > prod_max
                prod_max = prod_new; % store this as max
                optimal_A = A;
                optimal_B = B;
                %if counter > 10, break; end; % break after 20 iterations if even if the prod is rising
                
            elseif prod_new < prod_max
                
                counter = counter+1;     
            end
            
        end
        
        % breaking condition
        if  abs(prod_new - prod_max) <= .05 & counter > 3 
           if exist('optimal_A')==1 && exist('optimal_B')==1
            A = optimal_A;
            B = optimal_B; 
           end
           break
        else
            prod_old = prod_new;
            if counter > 3,
                if exist('optimal_A')==1 && exist('optimal_B')==1
                    A = optimal_A;
                    B = optimal_B; 
                end
                break; 
            end; % break after 20 iterations if even if the prod is rising
        end
    end
    level0.child0.prod(iter) = prod_new;                 
    iter = iter+1;                        
end % while 
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

%level(0).child(1)= struct('data',sD.data);6

%^^^^^^^^^^  BINARY TREE LEVEL 1 ONWOARDS ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


il =1; % initialize the level loop
while il <= 4 % loop till the breaking condition   
    try, % exception for the level error
        if isstruct(level(il))==1
            %fprintf ( 'Level \n',i)
            %clear children;
            cc = 1; % Initialize the child return
            for jl = 1:2^il
                try, 
                    %if isstruct(level(i).child(j))
                    if ~isempty(level(il).child(jl).data)%&(isstruct(level(i).child(j)))    
                        %s = length(level(i).child(j).data);
                        
                        fprintf(1,' level / child [%d, %d] \n',il, jl);   
                        children([cc]) = jl;
                        cc = cc+1;
                        
                        prod_old = 100*100/(100+100);
                        prod_max = 0;
                        counter = 0; counter0 = 0; max_counter = 0; iter = 1; counterNan =1;
                        
                        while 1
                            [centroids, clusters,err] = kmeans ('batch',level(il).child(jl).data,2, 'epochs',10000);
                            clusters = [(1:length(level(il).child(jl).data(:,1)))' clusters]; % annotate clusters with indix
                            
                            % Find out the indices of the vectors (rows) for both classes 
                            l= 0; k= 0;       
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
                            std_parent = mean(std(level(il).child(jl).data)); % take mean of parent (try different criterian, as well)
                            
                            % 1) 
                            %varDiff = abs(var (centroids(1,:)) - var(centroids(2,:))); % not good results
                            % comments ---- feels like it is not true representative of the children, stops too early
                            
                            % 2)
                            %varDiff = mean(var (level(il).child(jl).data(A,:))-var(level(i).child(jl).data(B,:)));
                            % comment ---- 
                            %----varDiff = abs(mean(var (level(i).child(jl).data(A,:))-var(level(i).child(jl).data(B,:))));
                            % comment ----- I think the best (on lung8)
                            
                            % 3)
                            %varDiff = mean(min(var (level(i).child(jl).data(A,:)),var(level(i).child(jl).data(B,:))));
                            %varDiff = mean(max(var (level(i).child(jl).data(A,:)),var(level(i).child(jl).data(B,:))));
                            % comment ---- MAX could not produce any result
                            %comments --- gave even bad results than 2
                            
                            % 4)
                            % This is the expresstion used in pat_cluster_intersection
                            % varDiff=  ( mean(varience(1,:))+  mean(varience(2,:)))/2
                            varDiff = (mean(var (level(il).child(jl).data(A,:)))+mean(var(level(il).child(jl).data(B,:))))/2;
                            %***************** END BREAKING CONDITION ********************************
                            
                            % break if the STD of parent becomes less that Difference of varience of the children.
                            if (std_parent < varDiff) ; level(il)>=4; break; end
                            %try
                               
                      
                                
                                
                                sD0 = som_data_struct (level(il).child(jl).data','comp_names',level(il).child(jl).label_names,'labels',sD.comp_names);                 
                                sD1 = som_data_struct (level(il).child(jl).data(A,:)','comp_names',level(il).child(jl).label_names(A),'labels',sD.comp_names);
                                sD2 = som_data_struct (level(il).child(jl).data(B,:)','comp_names',level(il).child(jl).label_names(B),'labels',sD.comp_names);
                                sM0 = som_make (sD0); 
                                %sM0 = checkTopol (sM0);
                                sM1 = som_make(sD1,'msize',toPology.topol.msize);%,'msize',sM.topol);
                                %sM1 = checkTopol (sM1);
                                sM2 = som_make(sD2,'msize',toPology.topol.msize);%,'msize',sM.topol);
                                %sM2 = checkTopol (sM2);
% Changes made in Nov 2003, to make run fast...this is origional                                
%                                 sM00 = som_autolabel_rank (sM0,sD0,'add1');
%                                 sM11 = som_autolabel_rank (sM1,sD1,'add1');
%                                 sM12 = som_autolabel_rank (sM2,sD2,'add1');         
                            try
                              % tic
                               %while toc < 100
                                sM00 = som_autolabel (sM0,sD0,'auto');
                                sM11 = som_autolabel (sM1,sD1,'auto');
                                sM12 = som_autolabel (sM2,sD2,'auto');      
                               %end
                            %catch
                             %   break
                            end
                            
                                % find the intersections and their indices in the parent codebook sM00 is the parent
                                [inds11,inds12]= label_comp(sD0.labels(:,1),sM11.labels(:,1),sM12.labels(:,1)); 
                                % find the probabilities of pridiction
                                [class_labels1,class_labels2] = gmm_bays_classify_3d(sD0,sM00,sD0.data(inds11,:),sD0.data(inds12,:));   
                                % the scores are the actual labels and the pridicted labels                            
                                score1 = 100*length(intersect(class_labels1,sD0.labels(inds11)))/length(inds11)
                                score2 = 100*length(intersect(class_labels2,sD0.labels(inds12)))/length(inds12)
                                
                                if score1 >= 95 & score2 >= 95; level(il).child(jl).prod(iter) = 3; level(il) >= 4; break;  end;
                                if isnan(score1)==1 | isnan(score2)==1, counterNan= counterNan+1
                                    if counterNan > 5
                                        if exist('optimal_A')==1 && exist('optimal_B')==1
                                            A = optimal_A;
                                            B = optimal_B; 
                                        end
                                        break; 
                                    end 
                                end
                                if score1==0 | score2==0, counter0= counter0+1
                                    if counter0 > 3
                                        if exist('optimal_A')==1 && exist('optimal_B')==1
                                            A = optimal_A;
                                            B = optimal_B; 
                                        end
                                        break; 
                                    end
                                end
                                if score1==0, score1=1; elseif score2==0, score2=1;end;
                                
                                % prod is the product of both scores divided by the sum
                                % Maximizing prod will garuntee optmal solution
                                if isnan(score1)==0 & isnan(score2)==0
                                    
                                    prod_new = (score1 * score2)/(score1 + score2) % maximum value of prod can take is 50
                                    
                                    
                                    if prod_new >= 40,level(il).child(jl).prod(iter) = prod_new; break; end;                                                                     
                                    % set the max_prod to maximum of prod_old and prod_new
                                    if prod_new <= 40 
                                        if prod_new > prod_max
                                            prod_max = prod_new; % store this as max
                                            optimal_A = A;
                                            optimal_B = B;
                                            %if counter > 10, break; end; % break after 20 iterations if even if the prod is rising
                                            
                                        elseif prod_new < prod_max
                                            
                                            counter = counter+1;     
                                        end
                                        
                                    end
                                end
                                % breaking condition
                                if  abs(prod_new - prod_max) <= .05 & counter > 3 
                                    if exist('optimal_A')==1 && exist('optimal_B')==1
                                        A = optimal_A;
                                        B = optimal_B; 
                                    end
                                    break
                                else
                                    prod_old = prod_new;
                                    if counter > 3,
                                        if exist('optimal_A')==1 && exist('optimal_B')==1
                                            A = optimal_A;
                                            B = optimal_B; 
                                        end
                                       
                                        break; 
                                    end; % break after 20 iterations if even if the prod is rising
                                end
                                %end
                            level(il).child(jl).prod(iter) = prod_new;                 
                            iter = iter+1;                        
                        end % 2n while
                        % record the prod for this level for plotting
                        
                        
                        level(il).child(jl).sD0 = sD0;
                        level(il).child(jl).sM0 = sM0;
                        level(il).child(jl).labelsA = A;
                        level(il).child(jl).labelsB = B;
                        level(il+1).child(2*jl-1).label_names = level(il).child(jl).label_names(A);
                        level(il+1).child(2*jl-1).data = level(il).child(jl).data(A,:);
                        level(il+1).child(2*jl-1).labels= level(il).child(jl).labels(A);  
                        level(il+1).child(2*jl-1).sD = sD1;
                        level(il+1).child(2*jl-1).sM = sM1;  
                        
                        level(il+1).child(2*jl).label_names=  level(il).child(jl).label_names(B);
                        level(il+1).child(2*jl).data = level(il).child(jl).data(B,:);
                        level(il+1).child(2*jl).labels= level(il).child(jl).labels(B);  
                        level(il+1).child(2*jl).sD= sD2;  
                        level(il+1).child(2*jl).sM= sM2;                       
                        
                        
                        if print == 1 
                            fid = fopen (filename,'a');
                            fprintf (fid,'\n%s\n   >> ',(['level(' int2str(il+1) ').child(' int2str(2*jl-1) ') ---> level(' ...
                                    int2str(il+2) ').child(' int2str(4*jl-3) ' & ' int2str(4*jl-2) ')']));
                            fprintf(fid, '\n');                               
                            for iii = A, fprintf (fid, '%s,   ',level(il).child(jl).label_names{iii});end 
                            fprintf(fid, '\n');
                            fclose (fid);
                            
                            fid = fopen (filename,'a');
                            fprintf (fid,'\n%s\n   >> ',(['level(' int2str(il+1) ').child(' int2str(2*jl) ') ---> level(' ...
                                    int2str(il+2) ').child(' int2str(4*jl-1) ' & ' int2str(4*jl) ')']));
                            fprintf(fid, '\n');
                            for iii = B, fprintf (fid, '%s,   ',level(il).child(jl).label_names{iii});end 
                            fprintf(fid, '\n');
                            fclose (fid);
                        end
                        
                        
                    end %if 
                    
                    
                catch, % if the i exceeds the generated level
                    %disp (lasterr) Never display this error as this error is wanted to go to next step
                    fff = 0; % simply take next j value
                    %break % break the for loop
                end % try
                
            end % for jl
            level(il).children = children;
            il =il+1; 
        else 
            
            break
        end % if
        
    catch, % if the i exceeds the generated level
        return % break the for loop and return the calling function
    end % try
end % for or while      