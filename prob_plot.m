
load probablities;
[prob_l, prob_r,prob_l_var,prob_r_var, prob_l_std,prob_r_std] = prob_gtree (gtree_classify)
clf
plot (prob_l, '-sb')
set(findobj(gca,'Type','line','color','b'),'color','b','LineWidth',3);  
hold on
plot (prob_r, '-dr')
set(findobj(gca,'Type','line','color','r'),'color','r','LineWidth',3);  



plot (prob_l-prob_l_var/2, '--b')
plot (prob_l+prob_l_var/2, '--b')
%set(findobj(gca,'Type','line','color','m'),'color','m','LineWidth',3);  
hold on
plot (prob_r-prob_r_var/2, '--r')
plot (prob_r+prob_r_var/2, '--r')
%set(findobj(gca,'Type','line','color','y'),'color','y','LineWidth',3);  

%plot (prob_l_std, '-sb')
%set(findobj(gca,'Type','line','color','m'),'color','m','LineWidth',3);  
%hold on
%plot (prob_r_std, '-dr')
%set(findobj(gca,'Type','line','color','y'),'color','y','LineWidth',3);  


title ('Average probablities of correct pridiction')
xlabel('Level of binary tree')
xlabel('Level of binary tree')
ylabel('probablity')
legend ('prob_l','-r')
legend ('Left Child','Right Child')
