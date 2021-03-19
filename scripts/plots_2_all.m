figure(2)
for(plot_number = 1 : howmany_systems)
    subplot(1,2,2)
    cc = jet(howmany_systems);
    
    ep_max_list_local = ep_max_list(:,plot_number);
    J_delta_ratio_local = J_delta_ratio(:,plot_number);
    
    plot(ep_max_list_local,J_delta_ratio_local,':ok','MarkerSize',15,'LineWidth',3,'color',cc(plot_number,:));
    hold on
    set(gca,'FontSize',25);
    set(gcf,'position',[10,10,1600,800])
    xlabel('$\epsilon$','Interpreter','latex','FontSize',25);
    ylabel('$({\hat{J}}^2-{J^*}^2)/{J^*}^2$','Interpreter','latex','FontSize',25);
    title('Suboptimality')
    xticks([0:0.2:2])
    legend({'$\rho$= 0.4','$\rho$= 0.5','$\rho$= 0.6','$\rho$= 0.7','$\rho$= 0.8','$\rho$= 0.9','$\rho$= 0.99'},'Interpreter','latex')
    axis([0 2 0 0.9])
    grid on
    
    subplot(1,2,1)
    plot(sigma_ok_list,ep_max_list_local,':ok','MarkerSize',15,'LineWidth',3,'color',cc(plot_number,:));
    hold on
    set(gca,'FontSize',25);
    set(gcf,'position',[10,10,1600,800])
    ylabel('$\epsilon$','Interpreter','latex','FontSize',25);
    xlabel('$\sigma$','Interpreter','latex','FontSize',25);
    title('Noise to estimation error')
    xticks([0:0.001:0.006])
    legend({'$\rho$= 0.4','$\rho$= 0.5','$\rho$= 0.6','$\rho$= 0.7','$\rho$= 0.8','$\rho$= 0.9','$\rho$= 0.99'},'Interpreter','latex')

    grid on
end