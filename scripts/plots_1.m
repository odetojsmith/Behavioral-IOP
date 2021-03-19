figure(1)
cc = jet(howmany_systems);
plot([0:len.n_horizon],y_real_data(:,1),'--o','MarkerSize',15,'LineWidth',3,'color',cc(plot_number,:)); 
hold on
scatter([0:len.n_horizon],y_real_data_lqg(:,1),100,'xk','LineWidth',3);

%figure(1)
%plot([0:len.n_horizon],y_real_data_KF(:,2),'-og','MarkerSize',13,'LineWidth',3);
%hold on
%plot([0:len.n_horizon],y_real_data(:,2),'--+k','MarkerSize',13,'LineWidth',3);
set(gca,'FontSize',30);
set(gcf,'position',[10,10,1000,750])
ylabel('System output $y_1$','Interpreter','latex','FontSize',30);
xlabel('$t$','Interpreter','latex','FontSize',30);
xticks([0:1:10])
grid on

h = legend('BIOP','LQG','Interpreter','latex');
set(h,'FontSize',30);

%axis([0 10 -0.01 0.5])
title('Optimal closed-loop trajectory')