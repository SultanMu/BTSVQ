function [child] = p_tree_btsvq (sD, sM, level,filename)

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

% Mujahid sultan, msultan@uhnres.utoronto.ca 
% beta 1.0
%############################################################################

  
  if ~isstruct(sM) 
    disp(['Map should be struct']);
    return
  end
% take the ranked labels from root level codebook
sM_rank = som_autolabel_rank (sM,sD,'add1');
%level(0).child(1).gene_labels = sM.labels;

i =1; % initialize the level loop
while 1 % loop till the breaking condition   
%for i = 1 : 8 % 50 levels of Partative tree (may be take it as argument at later stage)     
    try, % exception for the level error
      if isstruct(level(i))
            %fprintf ( 'Level \n',i)
            
            for j = 1:2:2^i
             try, 
                %if isstruct(level(i).child(j))
                if (level(i).child(j).data ~=[]) 
                    fprintf(1,' level / child [%d, %d] \n',i, j);   

                    %s = length(level(i).child(j).data);
                    labels_child = level(i).child(j).label_names; 
                    sD1 = som_data_struct (level(i).child(j).data','labels',sD.labels,'comp_names',level(i).child(j).label_names);
					sD2 = som_data_struct (level(i).child(j+1).data','labels',sD.labels,'comp_names',level(i).child(j+1).label_names);
					
                    sM11 = som_make_super(sD1);%,'msize',sM.topol);
                    sM11 = som_autolabel_rank (sM11,sD1,'add1');
                    sM12 = som_make_super(sD2);%,'msize',sM.topol);
                    sM12 = som_autolabel_rank (sM12,sD2,'add1');
                    

%                     clf % to avoid the screwed up figures
%                     som_show (sM11, 'compi', 'all','bar','none', 'footnote','' );
%                     % print the figure
%                     figno = char(strcat({'level_'}, {int2str(i)}, {'_child_'},{int2str(j)}))
%                     print (gcf, '-djpeg', '-r200', figno)
%                     clf
%                     som_show (sM12, 'compi', 'all','bar','none', 'footnote','' );
%                     % print the figure
%                     figno = char(strcat({'level_'}, {int2str(i)}, {'_child_'},{int2str(j+1)}))
%                     print (gcf, '-djpeg', '-r200', figno)
%                     
                    
					%    m11 = m0 intersection { c11 - (c11 intersection c12)}
					%    m12 = m0 intersection { c12 - (c11 intersection c12)}
					%    c11 - c12 = DD
                    
                    % Take the labels of parent for this child.
                    %if i == 1
                        m0 = sM_rank.labels(:,1);
                       [m11,inds11,m12,inds12]= gene_label_comp(m0,sM11.labels(:,1),sM12.labels(:,1)); % write this function
                    level(i).child(j).gene_labels = m11;
                    level(i).child(j+1).gene_labels = m12;                    
   
%                     else
%                         % Get the gene labels from the upper level and pass to the next
%                          %m0 = level(i-1).child((j+1)/2).gene_lables(:,1);                   
% 
%                          m0 = sM_rank.labels(:,1); % try this, as the upper did not gave any gene label after step 1
%                         [m11,inds11,m12,inds12]= gene_label_comp (m0,sM11.labels(:,1),sM12.labels(:,1)); % write this function
%                     
%                     level(i).child(j).gene_labels = m11;
%                     level(i).child(j+1).gene_labels = m12;                    
%                     end
                
                
                    clf; 
                    plot(sM.codebook(inds11,:)','-');                     
                    figno = char(strcat({'level_'}, {int2str(i)}, {'_child_'},{int2str(j)},{'_profile'}))
                    print (gcf, '-djpeg', '-r200', figno)
                    clf;
                    plot(sM.codebook(inds12,:)','-');
                    figno = char(strcat({'level_'}, {int2str(i)}, {'_child_'},{int2str(j+1)},{'_profile'}))
                    print (gcf, '-djpeg', '-r200', figno)
                    

                     fid = fopen (filename,'a');
                     fprintf (fid,'\n%s\n   >> ',(['level(' int2str(i) ').child(' int2str(j) ')'])); 
                     fprintf(fid, '\n');                               
                     for iii = 1:length(level(i).child(j).gene_labels)
                         fprintf (fid, '%s\t',char(level(i).child(j).gene_labels(iii)));
                     end 
                     fprintf(fid, '\n');
                     
                     
                     fprintf (fid,'\n%s\n   >> ',(['level(' int2str(i) ').child(' int2str(j+1) ')']));    
                     fprintf(fid, '\n');                               
                     for iii = 1:length (level(i).child(j+1).gene_labels)
                         fprintf (fid, '%s\t',char(level(i).child(j+1).gene_labels(iii)));
                     end 
                     fprintf(fid, '\n');
                     fclose (fid);

                    
                        %end
                end %if 
            catch, % if the i exceeds the generated level
                fff = 0; % simply take next j value
                %break % break the for loop
            end % try
        end % for j
        i =i+1; 
    else 
        break
    end % if
catch, % if the i exceeds the generated level
    return % break the for loop and return the calling function
end % try
end % for or while
