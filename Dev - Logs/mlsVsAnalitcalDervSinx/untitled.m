clc
close all

L= 3*2.032;
k = 2*pi/L;
udervA = [cos(k*x), -sin(k*x), -cos(k*x)];

figure(1)
hold on
for i=1:3    
    subplot('Position', [0.06,0.08+0.33*(3-i),0.93,0.21])
    hold on
    plot(x/L,udervA(:,i),'k','LineWidth',4)
    plot(x/L,uderv(:,i),'r-.','LineWidth',4)
    legend('Analytical','Numerical', 'location', 'eastoutside')    
    grid on
    set(gca,'GridAlpha',1,'GridLineStyle','--','FontSize',20)
    ylim([-1.5 1.5])
    xlim([-2 39]/L)
    
    xlabel('$x/L$', 'Interpreter', 'latex', 'FontSize', 20)
    
    if(i==1)
        title('1st Derivative of $\sin(2\pi\,x/L)$. Mesh resolution = $L/30$','Interpreter','Latex','FontSize',20)
        ylabel('$d / d x$','Interpreter','Latex','FontSize',20)
    elseif(i==2)
        title('2nd Derivative of $\sin(2\pi\,x/L)$. Mesh resolution = $L/30$','Interpreter','Latex','FontSize',20)
        ylabel('$d^2 / d x^2$','Interpreter','Latex','FontSize',20)
    elseif(i==3)
        title('3rd Derivative of $\sin(2\pi\,x/L)$. Mesh resolution = $L/30$','Interpreter','Latex','FontSize',20)
        ylabel('$d^3 / d x^3$','Interpreter','Latex','FontSize',20)
    end
        
end
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.05, 0.01, 0.9, 0.95]);