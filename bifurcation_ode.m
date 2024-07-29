clear;%close all;
global alpha delta mu phi A omega 

%% parameter(1):Vary with A
alpha=1.1;delta=0.8;mu=0.3;phi=0.15;omega=0.942;
% A=1.3;

%% parameter(2):Vary with phi
% alpha=0;delta=0.8;mu=0.1;A=1;omega=0.942;
% % phi=0.2;

%%
i=1;
AA=0.5:0.01:2;% 设置你需要计算的参数的范围 A
% AA=0:0.001:0.6;% 设置你需要计算的参数的范围 phi
x0=[displacement(1,1); displacement(2,1)];% 以TMRM的第一个解作为初值
% x0=[0.2; 0.1];

L=length(AA);
figure
for l=1:L
    A=AA(l);% 设置你需要计算的参数 A
    % phi=AA(l);% 设置你需要计算的参数 phi
    Tdata=0:0.01:500;
    options=odeset('RelTol',1e-10,'AbsTol',1e-10);
    [t,num] = ode45(@neuron,Tdata,x0,options);
    
    xmax_2=num(45001:50001,1);
    xmax(i).xmax=getmax(xmax_2);
    xmin(i).xmin=getmin(xmax_2);
    i=i+1;
    disp(i)
end
for j=1:i-1
    QQ=AA(j)*ones(1,length(xmax(j).xmax));
    p=plot(QQ,xmax(j).xmax,'k.','MarkerSize',5);
    hold on;
    % QQ1=AA(j)*ones(1,length(xmin(j).xmin));
    % plot(QQ1,xmin(j).xmin,'b.','MarkerSize',5);
    % hold on;
end
h1=legend(p,{'$$ Numerical $$'});
set(h1,'Interpreter','latex','FontSize',15);
% xlabel('$ A $','Interpreter','latex');
% ylabel('$ x_1 $','Interpreter','latex')
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);


%% 函数
function dydt = neuron(t,y)
global  alpha delta mu phi A omega
dydt=zeros(2,1);
dydt(1)=y(1)*(1-phi)-(1/3)*y(1).^3-y(2)+A*cos(omega*t);
dydt(2)=mu*(y(1)+alpha-delta*y(2));

end
