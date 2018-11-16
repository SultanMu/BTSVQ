function [] = plot_pridictions (ptree_classify,level0)

max_all = 0;
for i = 1:length(ptree_classify)
    maxi = length(ptree_classify(i).children)
    if maxi > max_all
        max_all = maxi
    end
end

figure;
subplot(length(ptree_classify),max_all,1)
%for l = 1:length(level0.child0.prod)
plot(level0.child0.prod,'-+b')
set(findobj(gca,'Type','line','Color','b'),'Color','blue','LineWidth',2);
%hold on            
%end

for i = 1:length(ptree_classify)-1
    
    for j = 1:length(ptree_classify(i).children)
        
        subplot(length(ptree_classify),max_all,(j+((i)*max_all)))
        set(gca,'XTickLabel',{''})
        if ~isempty(ptree_classify(i).child(ptree_classify(i).children(j)).prod)==1
            
%             for k = 1 :length(ptree_classify(i).child(ptree_classify(i).children(j)).prod)
%                 
%                 if isnan(ptree_classify(i).child(ptree_classify(i).children(j)).prod(k)) | (isempty(ptree_classify(i).child(ptree_classify(i).children(j)).prod(k)) ==1)
%                     ptree_classify(i).child(ptree_classify(i).children(j)).prod(k) = 1;
%                 end
%                 
%             end
            ptree_classify(i).child(ptree_classify(i).children(j)).prod
            plot(ptree_classify(i).child(ptree_classify(i).children(j)).prod,'-+b')
            set(findobj(gca,'Type','line','Color','b'),'Color','blue','LineWidth',2);
            hold on            
        end
        
        
    end
    
end

print (gcf, '-djpeg', '-r200', 'plot_of_pridictions');