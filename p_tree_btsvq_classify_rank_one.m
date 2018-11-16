function [child] = p_tree_btsvq_classify_rank_one (sD, sM, level,filename)

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
% activate this *****************************
[sM0,sM0_qerr]= som_autolabel_rank(sM,sD,'add1');
%sM_rank = som_autolabel_rank (sM,sD,'add1','ranked');
%sM_rank = som_autolabel(sM,sD,'vote');
%level(0).child(1).gene_labels = sM.labels;
scrsz = get(0,'ScreenSize');
figure(5);
set(5,'Position',[1/2 scrsz(4)/2 scrsz(3)/1.2 scrsz(4)/2]);

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
                        
                        
                        sD1 = level(i).child(j).sD;
                        sD2 = level(i).child(j+1).sD;
                        sM11 = som_autolabel_rank (level(i).child(j).sM,level(i).child(j).sD,'add1');
                        sM12 = som_autolabel_rank (level(i).child(j+1).sM,level(i).child(j+1).sD,'add1');                                        
                        
                        %    This is the fancy part                    
                        %    m11 = m0 intersection { c11 - (c11 intersection c12)}
                        %    m12 = m0 intersection { c12 - (c11 intersection c12)}
                        %    c11 - c12 = DD
                        
                        % Take the labels of parent for this child.
                        %if i == 1  % we can build a system for passing different codebook or the origional data to the children
                        if i == 1
                            labels1 = level(i).child(j).labels; 
                            labels2 = level(i).child(j+1).labels; 
                            m0 = sM0.labels(:,1);
                        else
                            sD0 = level(i-1).child(ceil(j/2)).sD0;
                            sM0 = level(i-1).child(ceil(j/2)).sM0;
                            sM0 = som_autolabel_rank(sM0,sD0, 'add1');
                            m0 = sM0.labels(:,1);
                            labels1 = level(i-1).child(ceil(j/2)).labelsA;
                            labels2 = level(i-1).child(ceil(j/2)).labelsB;
                        end
                        
                        %[m11,inds11,m12,inds12]= gene_label_comp_3(m0,sM11.labels(:,1),sM12.labels(:,1)); 
                        [m11,inds11,m12,inds12,c012,inds012]= gene_label_comp_3(m0,sM11.labels(:,1),sM12.labels(:,1)); 
                        % the following part is to check the labels in sD, it is ploted in blck thin line........make it off in final ver
                        [d11,indsd11,d12,indsd12,d012,indsd012]= gene_label_comp_3(sD.labels,sM11.labels(:,1),sM12.labels(:,1)); 
                        level(i).child(j).gene_labels = m11;
                        level(i).child(j+1).gene_labels = m12;
                        level(i).common_gene_labels = c012;
                        
                        genes11 = sM0.labels(inds11,1);
                        [dummy rank11] = sort (sM0_qerr.labels(inds11,1));
                        ranked_genes11 = genes11(rank11);
                        cated_ranked_genes11 = strcat ((ranked_genes11),'(' ,dummy,')');
                        
                        genes12 = sM0.labels(inds12,1);
                        [dummy rank12] = sort (sM0_qerr.labels(inds12,1));
                        ranked_genes12 = genes12(rank12);
                        cated_ranked_genes12 = strcat ((ranked_genes12),'(' ,dummy,')');
                        
                        if ~isempty (inds012) 
                            genes012 = sM0.labels(inds012,1);
                            [dummy rank012] = sort (sM0_qerr.labels(inds012,1));
                            ranked_genes012 = genes012(rank012);
                            cated_ranked_genes012 = strcat ((ranked_genes012),'(' ,dummy,')');
                        end
                        
                        
                                                
                        % the plot/printing section
                        % Note we can not plot the actual data as the inds11 and ands12 are the indices of the labels selected by the
                        % SOM's of children from the ranked labels of the root. But can be done by writing another function 
                        % for returning the indices of the children SOM labels from the root data set.
                        
                        % third plot option
                        % we will plot the genes selected by the SOM's for both children 
                        
                        %if handles.classification =='yes'
                        %                         [posterior, class] = gmm_bays_classify_3c(level(i).child(j).data',level(i).child(j+1).data',sD1,sM11, sD2,sM12);   
                        %                         class1 = (length(sM0.labels(posterior(:,1)>.8,1))/length (sM0.labels))*100
                        %                         class2 = (length(sM0.labels(posterior(:,2)>.8,2))/length (sM0.labels))*100
                        
                        %end
                        clf; 
                        if inds11 ~= [] & length(inds11)>1
                            figure(5);
                            clf;
                            plot(1:length(labels1),(sM0.codebook(inds11,labels1)),'-c',...
                                1:length(labels1),mean(sM0.codebook(inds11,labels1)),'-ob', ...
                                (length(labels1)+1):(length(labels1)+length(labels2)),(sM0.codebook(inds11,labels2)),'-g', ...
                                (length(labels1)+1):(length(labels1)+length(labels2)),mean(sM0.codebook(inds11,labels2)),'-ob');
                            set(findobj(gca,'Type','line','Color','b'),'Color','blue','LineWidth',3);
                            hold on
                        
                            if ~isempty (indsd012)
                                plot (1:length(labels1),sD.data(indsd012,labels1),'-r',...
                                    (length(labels1)+1):(length(labels1)+length(labels2)),sD.data(indsd012,labels2),'-r');
                                %set(findobj(gca,'Type','line','Color','r'),'Color','r','LineWidth',2);
                                    %legend (d012,-1);    % make it off
                            end
                            hold on;
                            %set(findobj(gca,'Type','line','Color',[0 0 1]),'Color','blue','LineWidth',3);
                            
                            %                         hold on;
                            %                         
                            %                         % find first ranked gene from ranked_genes11  in sM0.codebook and plot its profile
                            %                         % May be add some dialog
                            %                         plot (1:length (labels1),sM0.codebook(strcmp(sM0.labels(:,1),ranked_genes11(1)),labels1),'-m', ...
                            %                             (length(labels1)+1):(length(labels1)+length(labels2)),sM0.codebook(strcmp(sM0.labels(:,1),ranked_genes12(1)),labels2),'-m');                               
                            %                         set(findobj(gca,'Type','line','Color','m'),'Color','m','LineWidth',2);
                            %                         % legend
                            %                         if length(cated_ranked_genes11) > 21
                            %                             legend (cated_ranked_genes11(1:21,1),-1);
                            %                         else
                            %                             legend (cated_ranked_genes11,-1);
                            %                         end
                            %                         if ~isempty (inds012)
                            %                             plot (1:length(labels1),sM0.codebook(inds012,labels1),'-k',...
                            %                                 (length(labels1)+1):(length(labels1)+length(labels2)),sM0.codebook(inds012,labels2),'-k');
                            %                             legend (cated_ranked_genes012,-1);    % make it off
                            %                         end
                            % print to the file
                                                      
                            clear tt;
                            clear pp;
                            for ij = 1:length (inds11)
                                t =  (mean (sM0.codebook(inds11(ij),labels2)) - mean(sM0.codebook(inds11(ij),labels1)));
                                de = (((var(sM0.codebook(inds11(ij),labels1)))^2 ./ length (labels1)) +((var(sM0.codebook(inds11(ij),labels2)))^2 ./ length (labels2)))^.5;
                                tt(ij) = t./de;
                                 [h,p,ci,stats] = ttest2(sM0.codebook(inds11(ij),labels2),sM0.codebook(inds11(ij),labels1));
                                 pp(ij) = p;
                            end
                            [max ranki] = sort(tt);
                            plot(1:length(labels1),(sM0.codebook(inds11(ranki(length(ranki))),labels1)),'-m',... 
                                    (length(labels1)+1):(length(labels1)+length(labels2)),(sM0.codebook(inds11(ranki(length (ranki))),labels2)),'-m');
                                set(findobj(gca,'Type','line','color','m'),'color','m','LineWidth',5);    
                            legend (sM0.labels(inds11(ranki(length(ranki)))),-1);
                            hold off;
                            % bring the variable deciding the specimen or Genes and then print the labels
                            xlabel('Specimens','fontsize',14);
                            ylabel('Normalized Expression Level');
                            title(char(strcat({'level\_'}, {int2str(i)}, {'   child\_'},{int2str(j)})),'Color','b','FontWeight','bold') 
                            figno = char(strcat({'level_'}, {int2str(i)}, {'_child_'},{int2str(j)},{'_profile'}))
                            print (gcf, '-djpeg', '-r200', figno)
                           
                            % Plot the t-stastisic
                            figure(6);
                            clf
                            subplot(3,2,3), plot (tt, '-or');
                            title (char(strcat({'level'}, {int2str(i)}, {'child'},{int2str(j)},{'(t-test)'})));
                            subplot(3,2,5), plot (pp, '-or');
                            title (char(strcat({'level'}, {int2str(i)}, {'child'},{int2str(j)},{'(p-test)'})));
                            %figno = char(strcat({'level_'}, {int2str(i)}, {' child_'},{int2str(j)},{'(t_test_p)'}));
                            %print (gcf, '-djpeg', '-r200', figno);
                            
                            % plot of probablitite 
                            % 1) take two children at any node and learn the mixtrue model form SM, and sD.
                            %    Classify all genes (taking the sub-sets of patients, to keep the dimentionality same)
                            %[posterior, class] = gmm_bays_classify_3c(level(i).child(j).data',level(i).child(j+1).data',sD1,sM11, sD2,sM12);   
                            %
                            % 2) take both children at a node and learn mixture model from sM and sD
                            %    Classify the root codebook (taking the sub-stes of the patients)
                            [posterior, class1_labels,class2_labels] = gmm_bays_classify_3c(sM0.labels(:,1),sM0.codebook(:,labels1),sM0.codebook(:,labels2),sD1,sM11, sD2,sM12);   
                            %figure(7)
                            %clf
                            subplot (3,2,1),plot (posterior(:,1),'.b');
                            hold on;
                            plot (posterior(:,2),'.r');
                            title (char(strcat({'level'}, {int2str(i)}, {'child'},{int2str(j)},{'(Prob)'})));
                            hold off;
                            subplot (3,2,2),plot (posterior(:,1),'.b');
                            hold on;
                            plot (posterior(:,2),'.r');
                            title (char(strcat({'level'}, {int2str(i)}, {'child'},{int2str(j+1)},{'(Prob)'})));
                            hold off;
                            hold on; %for figure 6
                            %figno = char(strcat({'level_'}, {int2str(i)}, {' child_'},{int2str(j)},{' & '},{int2str(j+1)},{'(prob)'}));
                            %print (gcf, '-djpeg', '-r200', figno);
                            
                            % class1_labels are the labels of sM0.codebook which were generated by mixture model_1 
                            % class2_labels are the labels of sM0.codebook which were generated by mixture model_2 
                            % Now the task is to check how much percent of class1_labels are in sM1.labels and Class2_labels are in sM2.labels
                            try 
                                [class1_in_sM0] = gene_label_comp_4 (sM0.labels(:,1),class1_labels);
                                [common1] = gene_label_comp_4 (class1_in_sM0,ranked_genes11);
                                level(i).child(j).percent  = (length(common1)/length (ranked_genes11))*100;
                                
                                % similarly
                                [class2_in_sM0] = gene_label_comp_4 (sM0.labels(:,1),class2_labels);
                                [common2] = gene_label_comp_4 (class2_in_sM0,ranked_genes11);
                                level(i).child(j+1).percent = (length(common2)/length(ranked_genes12))*100;
                            catch
                                level(i).child(j).percent = 0;
                                level(i).child(j+1).percent =0;
                                
                            end
                            
                            
                        elseif length(inds11) == 1
                            figure(5);
                            clf;
                            plot(1:length(labels1),(sM0.codebook(inds11,labels1)),'-c',...
                                (length(labels1)+1):(length(labels1)+length(labels2)),(sM0.codebook(inds11,labels2)),'-g')
                                                      
                            set(findobj(gca,'Type','line','Color','c'),'Color','blue','LineWidth',3);
                            %set(findobj(gca,'Type','line','Color',[0 0 1]),'Color','blue','LineWidth',3);
                            title(char(strcat({'level\_'}, {int2str(i)}, {'   child\_'},{int2str(j)})),'Color','b','FontWeight','bold') 
                            % legend
                            if length(cated_ranked_genes11) > 21
                                legend (cated_ranked_genes11(1:21,1),-1);
                            else
                                legend (cated_ranked_genes11,-1);
                            end
                            if ~isempty (inds012)
                                plot (1:length(labels1),sM0.codebook(inds012,labels1),'-k',...
                                    (length(labels1)+1):(length(labels1)+length(labels2)),sM0.codebook(inds012,labels2),'-k');
                                    set(findobj(gca,'Type','line','Color','k'),'Color','k','LineWidth',2);
                                legend (cated_ranked_genes012,-1);    % make it off
                            end
                            
                            % print to the file
                            figno = char(strcat({'level_'}, {int2str(i)}, {'_child_'},{int2str(j)},{'_profile'}))
                            print (gcf, '-djpeg', '-r200', figno)
                            
                            
                            % plot of probablitite 
                            % 1) take two children at any node and learn the mixtrue model form SM, and sD.
                            %    Classify all genes (taking the sub-sets of patients, to keep the dimentionality same)
                            %[posterior, class] = gmm_bays_classify_3c(level(i).child(j).data',level(i).child(j+1).data',sD1,sM11, sD2,sM12);   
                            %
                            % 2) take both children at a node and learn mixture model from sM and sD
                            %    Classify the root codebook (taking the sub-stes of the patients)
                            
                            % plot of probablitite 
                            % 1) take two children at any node and learn the mixtrue model form SM, and sD.
                            %    Classify all genes (taking the sub-sets of patients, to keep the dimentionality same)
                            %[posterior, class] = gmm_bays_classify_3c(level(i).child(j).data',level(i).child(j+1).data',sD1,sM11, sD2,sM12);   
                            %
                            % 2) take both children at a node and learn mixture model from sM and sD
                            %    Classify the root codebook (taking the sub-stes of the patients)
                            [posterior, class1_labels,class2_labels] = gmm_bays_classify_3c(sM0.labels(:,1),sM0.codebook(:,labels1),sM0.codebook(:,labels2),sD1,sM11, sD2,sM12);   
                            figure(7);
                            clf;
                            plot (posterior(:,1),'.b')
                            hold on
                            plot (posterior(:,2),'.r')
                            figno = char(strcat({'level_'}, {int2str(i)}, {' child_'},{int2str(j)},{' & '},{int2str(j+1)},{'(prob)'}))
                            print (gcf, '-djpeg', '-r200', figno)
                            
                            % class1_labels are the labels of sM0.codebook which were generated by mixture model_1 
                            % class2_labels are the labels of sM0.codebook which were generated by mixture model_2 
                            % Now the task is to check how much percent of class1_labels are in sM1.labels and Class2_labels are in sM2.labels
                            try
                                [class1_in_sM0] = gene_label_comp_4 (sMo.labels(:,1),class1_labels);
                                [common1] = gene_label_comp_4 (class1_in_sM0,ranked_genes11);
                                level(i).child(j).percent  = (length(common1)/length (ranked_genes11))*100;
                                
                                % similarly
                                [class2_in_sM0] = gene_label_comp_4 (sMo.labels(:,1),class2_labels);
                                [common2] = gene_label_comp_4 (class2_in_sM0,ranked_genes11);
                                level(i).child(j+1).percent = (length(common2)/length(ranked_genes12))*100;
                            catch
                                level(i).child(j).percent = 0;
                                level(i).child(j+1).percent =0;
                            end
                        end
                        
        % other child  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        if inds12 ~= [] & length(inds12)>1
                        
                            
                            figure(5);
                            clf;
                            plot(1:length(labels1),(sM0.codebook(inds12,labels1)),'-g',... 
                                1:length(labels1),mean(sM0.codebook(inds12,labels1)),'-ob', ...
                                (length(labels1)+1):(length(labels1)+length(labels2)),(sM0.codebook(inds12,labels2)),'-c',...
                                (length(labels1)+1):(length(labels1)+length(labels2)),mean(sM0.codebook(inds12,labels2)),'-ob');
                            set(findobj(gca,'Type','line','Color','b'),...
                                'Color','blue',...
                                'LineWidth',3);
                            hold on
                            if ~isempty (indsd012)
                                plot (1:length(labels1),sD.data(indsd012,labels1),'-r',...
                                    (length(labels1)+1):(length(labels1)+length(labels2)),sD.data(indsd012,labels2),'-r');
                                %set(findobj(gca,'Type','line','Color','r'),'Color','r','LineWidth',2);
                                    %legend (d012,-1);    % make it off
                            end
                            
                            %                         hold on;
                            %                         
                            %                         % find first ranked gene from ranked_genes11  in sM0.codebook and plot its profile
                            %                         % May be add some dialog
                            %                         plot (1:length (labels1),sM0.codebook(strcmp(sM0.labels(:,1),ranked_genes12(1)),labels1),'-m', ...
                            %                             (length(labels1)+1):(length(labels1)+length(labels2)),sM0.codebook(strcmp(sM0.labels(:,1),ranked_genes11(1)),labels2),'-m');                               
                            %                         set(findobj(gca,'Type','line','Color','m'),'Color','m','LineWidth',2);
                            %                         % legend
                            %                         if length(cated_ranked_genes12) > 21
                            %                             legend (cated_ranked_genes12(1:21,1),-1);
                            %                         else
                            %                             legend (cated_ranked_genes12,-1);
                            %                         end
                            %                         
                            %                         if ~isempty (inds012)
                            %                             plot (1:length(labels1),sM0.codebook(inds012,labels1),'-k',...
                            %                                 (length(labels1)+1):(length(labels1)+length(labels2)),sM0.codebook(inds012,labels2),'-k');
                            %                             legend (cated_ranked_genes012,-1);    % make it off
                            %                         end
                            try  % put this in try catch as some times the tt may get too smaller and we may get the out of dim error
                                clear tt;
                                clear pp;
                                for ij = 1:length (inds12)
                                    t =  (mean (sM0.codebook(inds12(ij),labels2)) - mean(sM0.codebook(inds12(ij),labels1)));
                                    de = (((var(sM0.codebook(inds12(ij),labels1)))^2 ./ length (labels1)) +((var(sM0.codebook(inds12(ij),labels2)))^2 ./ length (labels2)))^.5;
                                    tt(ij) = t./de;
                                    [h,p,ci,stats] = ttest2(sM0.codebook(inds12(ij),labels1),sM0.codebook(inds12(ij),labels2));
                                    pp(ij) = p;
                                end
                                
                                [max ranki] = sort(tt);
                                plot(1:length(labels1),(sM0.codebook(inds12(ranki(length(ranki))),labels1)),'-m',... 
                                    (length(labels1)+1):(length(labels1)+length(labels2)),(sM0.codebook(inds12(ranki(length (ranki))),labels2)),'-m');
                                set(findobj(gca,'Type','line','color','m'),'color','m','LineWidth',5);    
                                legend (sM0.labels(inds11(ranki(length(ranki)))),-1);
                            catch
                                disp (lasterr);
                            end
                            xlabel('Specimens','fontsize',14);
                            ylabel('Normalized Expression Level');
                            title(char(strcat({'level\_'}, {int2str(i)}, {'   child\_'},{int2str(j+1)})),'Color','b','FontWeight','bold')                             
                            hold off;
                            figno = char(strcat({'level_'}, {int2str(i)}, {'_child_'},{int2str(j+1)},{'_profile'}))
                            print (gcf, '-djpeg', '-r200', figno)
                            
                                                        
                            % plot of probablitite 
                            % 1) take two children at any node and learn the mixtrue model form SM, and sD.
                            %    Classify all genes (taking the sub-sets of patients, to keep the dimentionality same)
                            %[posterior, class] = gmm_bays_classify_3c(level(i).child(j).data',level(i).child(j+1).data',sD1,sM11, sD2,sM12);   
                            %
                            % 2) take both children at a node and learn mixture model from sM and sD
                            %    Classify the root codebook (taking the sub-stes of the patients)
                            % plot of probablitite 
                            % 1) take two children at any node and learn the mixtrue model form SM, and sD.
                            %    Classify all genes (taking the sub-sets of patients, to keep the dimentionality same)
                            %[posterior, class] = gmm_bays_classify_3c(level(i).child(j).data',level(i).child(j+1).data',sD1,sM11, sD2,sM12);   
                            %
                            % 2) take both children at a node and learn mixture model from sM and sD
                            %    Classify the root codebook (taking the sub-stes of the patients)
                            
                            % class1_labels are the labels of sM0.codebook which were generated by mixture model_1 
                            % class2_labels are the labels of sM0.codebook which were generated by mixture model_2 
                            % Now the task is to check how much percent of class1_labels are in sM1.labels and Class2_labels are in sM2.labels

                            
                            figure(6);
                            subplot(3,2,4), plot (tt, '-ob');
                            title (char(strcat({'level'}, {int2str(i)}, {'child'},{int2str(j+1)},{'(t-test)'})));
                            subplot(3,2,6), plot (pp, '-ob');
                            title (char(strcat({'level'}, {int2str(i)}, {'child'},{int2str(j+1)},{'(p-test)'})));
                            hold off;                            
                            figno = char(strcat({'level_'}, {int2str(i)}, {' child_'},{int2str(j)},{'& child_'},{int2str(j+1)},{'(statistic)'}));
                            print (gcf, '-djpeg', '-r200', figno);
                            
                        elseif length(inds12) == 1
                            figure(5);
                            clf;
                            plot(1:length(labels1),(sM0.codebook(inds12,labels1)),'-g',... 
                                (length(labels1)+1):(length(labels1)+length(labels2)),(sM0.codebook(inds12,labels2)),'-c');
                            
                            set(findobj(gca,'Type','line','Color','g'),...
                                'Color','blue',...
                                'LineWidth',3);
                            
                            title(char(strcat({'level\_'}, {int2str(i)}, {'   child\_'},{int2str(j+1)})),'Color','b','FontWeight','bold') 
                            %                         if length(cated_ranked_genes12) > 21
                            %                             legend (cated_ranked_genes12(1:21,1),-1);
                            %                         else
                            %                             legend (cated_ranked_genes12,-1);
                            %                         end
                            %                         if ~isempty (inds012)
                            %                             plot (1:length(labels1),sM0.codebook(inds012,labels1),'-k',...
                            %                                 (length(labels1)+1):(length(labels1)+length(labels2)),sM0.codebook(inds012,labels2),'-k');
                            %                             legend (cated_ranked_genes012,-1);    % make it off
                            %                         end
                            
                            figno = char(strcat({'level_'}, {int2str(i)}, {' child_'},{int2str(j)},{' & '},{int2str(j+1)},{'(prob)'}));
                            print (gcf, '-djpeg', '-r200', figno)
                            
                            % plot of probablitite 
                            % 1) take two children at any node and learn the mixtrue model form SM, and sD.
                            %    Classify all genes (taking the sub-sets of patients, to keep the dimentionality same)
                            %[posterior, class] = gmm_bays_classify_3c(level(i).child(j).data',level(i).child(j+1).data',sD1,sM11, sD2,sM12);   
                            %
                            % 2) take both children at a node and learn mixture model from sM and sD
                            %    Classify the root codebook (taking the sub-stes of the patients)
                            % plot of probablitite 
                            % 1) take two children at any node and learn the mixtrue model form SM, and sD.
                            %    Classify all genes (taking the sub-sets of patients, to keep the dimentionality same)
                            %[posterior, class] = gmm_bays_classify_3c(level(i).child(j).data',level(i).child(j+1).data',sD1,sM11, sD2,sM12);   
                            %
                            % 2) take both children at a node and learn mixture model from sM and sD
                            %    Classify the root codebook (taking the sub-stes of the patients)
                            
                        end
                        
                        % add the cated gene labels to the structure
                        level(i).child(j).ranked_gene_labels = cated_ranked_genes11;
                        level(i).child(j+1).ranked_gene_labels =cated_ranked_genes12;                    
                        
                        try % avoid any jump, by avoiding mistake in writing, may be empty string                  
                            % the g-tree text file section
                            fid = fopen (filename,'a');
                            fprintf(fid, ' Correct Percentage = %s\n', int2str(level(i).child(j).percent));
                            fprintf (fid,'\n%s\n   >> ',(['level(' int2str(i) ').child(' int2str(j) ')'])); 
                            fprintf(fid, '\n');                               
                            for iii = 1:length(level(i).child(j).ranked_gene_labels)
                                fprintf (fid, '%s\t',char(level(i).child(j).ranked_gene_labels(iii)));
                            end 
                            
                            fprintf(fid, '\n');
                            fprintf(fid, ' Correct Percentage = %s\n', int2str(level(i).child(j+1).percent));
                            fprintf (fid,'\n%s\n   >> ',(['level(' int2str(i) ').child(' int2str(j+1) ')']));    
                            fprintf(fid, '\n');                               
                            for iii = 1:length (level(i).child(j+1).ranked_gene_labels)
                                fprintf (fid, '%s\t',char(level(i).child(j+1).ranked_gene_labels(iii)));
                            end 
                            fprintf(fid, '\n');
                            fclose (fid);
                        catch
                            fclose(fid); % close any axcedental opend file
                        end
                        
                        %end % if 
                    end %if 
                catch, % if the i exceeds the generated level
                    disp(lasterr)
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
