figure(2)
subplot(1,2,2)
cc = jet(howmany_systems);

ep_max_list_local = ep_max_list(:,plot_number);
J_delta_ratio_local = J_delta_ratio(:,plot_number);

plot(ep_max_list_local,J_delta_ratio_local,':ok','MarkerSize',10,'LineWidth',3,'color',cc(plot_number,:));
hold on
set(gca,'FontSize',18);
set(gcf,'position',[10,10,600,400])
xlabel('$\epsilon$','Interpreter','latex','FontSize',25);
ylabel('$({\hat{J}}^2-{J^*}^2)/{J^*}^2$','Interpreter','latex','FontSize',18);
grid on

subplot(1,2,1)
plot(sigma_ok_list,ep_max_list_local,':ok','MarkerSize',10,'LineWidth',3,'color',cc(plot_number,:));
hold on
set(gca,'FontSize',18);
set(gcf,'position',[10,10,600,400])
ylabel('$\epsilon$','Interpreter','latex','FontSize',30);
xlabel('$\sigma$','Interpreter','latex','FontSize',25);
grid on
