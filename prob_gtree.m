function [prob_l_mean, prob_r_mean,prob_l_var,prob_r_var,prob_l_std,prob_r_std] = prob_gtree (level)

% Partitions patient space using partitive clustering
% Generates SOM of each child
% Takes ranked gene labels form the respective SOM
% Computes m11 and m12 from the m0 from upper step
% Passes m0 to all children to compute m11 and m12 for that child
% plots profile of m11 and m12 but from only the codebook, mkae it to plot for the origional data.


%
%  INPUTS
% sM = data struct with fields
%       sM.data
%       sM.comp_names  
%       sM.labels;  labels of the data vectors
%
% ptree = data struct generated with  function (p_tree)
%
% OUTPUTS
% child = data struct with fild child containing labels for hthat child
% Child(i).labels

% written by Mujahid sultan, msultan@uhnres.utoronto.ca 
% beta 1.2, Sep 2002
%############################################################################




    prob_l(1) = max(level(1).child(1).prob, level(1).child(2).prob);
    prob_r(1) = max(level(1).child(1).prob, level(1).child(2).prob);

    prob_l_mean(1)= mean([level(1).child(1).prob level(1).child(2).prob]);
    prob_r_mean(1)= mean([level(1).child(1).prob level(1).child(2).prob]);
    prob_l_var(1)= var([level(1).child(1).prob level(1).child(2).prob]);
    prob_r_var(1)= var([level(1).child(1).prob level(1).child(2).prob]);
    prob_l_std(1)= std([level(1).child(1).prob level(1).child(2).prob]);
    prob_r_std(1)= std([level(1).child(1).prob level(1).child(2).prob]);
    k=1;
    l=1;
    i =2; % initialize the level loop
    while 1 % loop till the breaking condition   
        %for i = 1 : 8 % 50 levels of Partative tree (may be take it as argument at later stage)     
        try, % exception for the level error
            if isstruct(level(i))
                %fprintf ( 'Level \n',i)
                
                for j = 1:2:2^i
                    try, 
                        %if isstruct(level(i).child(j))
                        if j< 2^i/2 
                            prob_left(k) = max(level(i).child(j).prob, level(i).child(j+1).prob);
                            k = k+1;
                        else
                            prob_right(l) = max(level(i).child(j).prob, level(i).child(j+1).prob);
                            l = l+1;
                        end %if 
                        

                    catch, % if the i exceeds the generated level
                        % disp(lasterr) *********** Never display this error this is just dummy
                        fff = 0; % simply take next j value
                        %break % break the for loop
                    end % try
                end % for j
                
%                 if length(prob_left)>1
%                     boxplot (prob_left)
%                     set(findobj(gca,'Type','line','color','b'),'color','b','LineWidth',3);  
%                 end
%                 hold on
%                 if length(prob_right)>1
%                     boxplot (prob_right)
%                 end
                %set(findobj(gca,'Type','line','color','r'),'color','r','LineWidth',3);  
                %prob_l(i) = (prob_left);
                %prob_r(i)=  (prob_right);
                prob_l_mean(i) = mean (prob_left);
                prob_r_mean(i)= mean (prob_right);
                
                prob_l_var(i) = var (prob_left);
                prob_r_var(i)= var (prob_right);

                prob_l_std(i) = std (prob_left);
                prob_r_std(i)=  std (prob_right);

                
                i =i+1; 
                
            else 
                break
            end % if
        catch, % if the i exceeds the generated level
            return % break the for loop and return the calling function
        end % try
        
    end % for or while
    
