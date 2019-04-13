function [level] = p_tree (sD,labels, filename,print)

% level  = p_tree(sDt, 'labels', 'filename.txt',0 or 1)
% level = data structure containing k-means results
% ith level contains 2^i children
% labels field is empty string (for future use)
% Level(1).child(1).data 
% Level(2).child(1).data
% 
% sD, sM = data struct and trained map struct
% dim = number of columens in the data set
% dlenght = number of rows in the data set
% Also See
% p_tree_parent2
% child = p_treePparent2 (sM,ptree)

% Mujahid sultan, msultan@uhnres.utoronto.ca 
% beta 1.0
%


% if isstruct(sD) 
%  data = sD.name; 
%  labels = sD.labels;
%  comp_names = D.comp_names;
%else 
%  data = inputname(1);
%  labels = inputname(2);
%end

% Divide the data set into two children using k-means
% [centroids,clusters,err] = kmeans(method, D, n, [epochs], [verbose])
% centroids are in the 'dim' dimention
[centroids, clusters,err] = kmeans ('batch',sD,2,'epochs',10000);
%[centroids, clusters,err] = kmeans ('',sD,2,'epochs',10000);
clusters = [(1:length(sD.data(:,1)))' clusters]; % annotate clusters with indix

     
     % Find out the indices of the vectors (rows) for both classes 
     % Note that you can do it quite easily but this is the initial implementation, so I keep it
     l= 0;
     k= 0;
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
     level(1).child(1).data= sD.data(A,:);
     level(1).child(1).labels= [A];%sD.labels(A); % store vector(rows) labels as well
     level(1).child(1).label_names= sD.labels(A);
          
     level(1).child(2).data = sD.data(B,:);
     level(1).child(2).labels= [B];%sD.labels(A); % store vector(rows) labels as well
     level(1).child(2).label_names= sD.labels(B);     
     
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
        if isstruct(level(i))==1
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
                if ~isempty(level(i).child(j).data) %&(isstruct(level(i).child(j)))    
                    %s = length(level(i).child(j).data);
                
                    fprintf(1,' level / child [%d, %d] \n',i, j);   
                    children([cc]) = j;
                    cc = cc+1;
                    
                    [centroids, clusters,err] = kmeans ('batch',level(i).child(j).data,2);
			        clusters = [(1:length(level(i).child(j).data(:,1)))' clusters]; % annotate clusters with indix
          			clear A;
                    clear B;
         			% Find out the indices of the vectors (rows) for both classes 
           	    	l= 0;
           		    k= 0;
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
                    if (std_parent > varDiff) 
                   
                                  
                       level(i+1).child(2*j-1).label_names = level(i).child(j).label_names(A);
                       level(i+1).child(2*j-1).data = level(i).child(j).data(A,:);
           	           level(i+1).child(2*j-1).labels= level(i).child(j).labels(A);  
                                  	           
                       level(i+1).child(2*j).label_names=  level(i).child(j).label_names(B);
                       level(i+1).child(2*j).data = level(i).child(j).data(B,:);
                       level(i+1).child(2*j).labels= level(i).child(j).labels(B);  
                    
                       
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
                 
                    else
                     
                     level(i+1).child(2*j-1).data = [];
                     level(i+1).child(2*j).data =   [];
                        %Level(i+1).child(2*j-1).lables = [];
                        %Level(i+1).child(2*j).labels = [];    
                    end
            
                else
                    % go and try next j, if all j th children are empty then exit
                    fff= 0; 
                    %level(i+1).child(2*j-1).data = [];
   	                %level(i+1).child(2*j).data =   [];
                end %if 
 
            catch, % if the i exceeds the generated level
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