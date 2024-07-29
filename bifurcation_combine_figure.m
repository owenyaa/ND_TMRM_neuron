clear;
% % close all;

figure
load 'bifurcation_unstable_period2_N=50_omega_new=0.471_A_1_0.01_2.mat'
N_dof=2;N_harm=50;
w=0.942;w_new=w/2;

for ii=1:1:47 % 1——1.47
    parameter_a=every_a(ii).parameter_a;
    A(ii)=every_a(ii).A;
    Tdata1=0:0.01:500;
    % w0=parameter_a(1,1);
    Harm_parameter_a=parameter_a(2:end,:);
    %% 计算方程残差
    x=zeros(N_dof,length(Tdata1));
    dx=zeros(N_dof,length(Tdata1));
    ddx=zeros(N_dof,length(Tdata1));
    % 有X,Y两个自由度
    % for i=1:N_harm   % i=1,3,5
    for j=1:N_dof
        for i=1:N_harm   % i=1,3,5
            x(j,:)=x(j,:)+Harm_parameter_a(i,2*j-1)*cos(i*w_new*Tdata1)+Harm_parameter_a(i,2*j)*sin(i*w_new*Tdata1);
%             dx(j,:)=dx(j,:)-w0*(2*i-1)*Harm_parameter_a(i,2*j-1)*sin((2*i-1)*w0*Tdata)+w0*(2*i-1)*Harm_parameter_a(i,2*j)*cos((2*i-1)*w0*Tdata);
%             ddx(j,:)=ddx(j,:)-(w0*(2*i-1))^2*Harm_parameter_a(i,2*j-1)*cos((2*i-1)*w0*Tdata)-(w0*(2*i-1))^2*Harm_parameter_a(i,2*j)*sin((2*i-1)*w0*Tdata);
        end
        x(j,:)=x(j,:)+parameter_a(1,2*j-1);
    end
    x=x';
    xmax_2=x(30001:50001,1);
    xmax(ii).xmax=getmax(xmax_2);
    xmin(ii).xmin=getmin(xmax_2);
    hold on;
    QQ=A(ii)*ones(1,length(xmax(ii).xmax));
    p1=plot(QQ,xmax(ii).xmax,'b.','MarkerSize',4);
    hold on;
end

hold on
%%
load 'bifurcation_unstable_period1_N=50_omega_new=0.471_A_0.93_0.01_1.47.mat'
for ii=2:1:54 % 0.94——1.47
    parameter_a=every_a(ii).parameter_a;
    A(ii)=every_a(ii).A;
    Tdata1=0:0.01:500;
    % w0=parameter_a(1,1);
    Harm_parameter_a=parameter_a(2:end,:);
    %% 计算方程残差
    x=zeros(N_dof,length(Tdata1));
    dx=zeros(N_dof,length(Tdata1));
    ddx=zeros(N_dof,length(Tdata1));
    % 有X,Y两个自由度
    % for i=1:N_harm   % i=1,3,5
    for j=1:N_dof
        for i=1:N_harm   % i=1,3,5
            x(j,:)=x(j,:)+Harm_parameter_a(i,2*j-1)*cos(i*w_new*Tdata1)+Harm_parameter_a(i,2*j)*sin(i*w_new*Tdata1);
%             dx(j,:)=dx(j,:)-w0*(2*i-1)*Harm_parameter_a(i,2*j-1)*sin((2*i-1)*w0*Tdata)+w0*(2*i-1)*Harm_parameter_a(i,2*j)*cos((2*i-1)*w0*Tdata);
%             ddx(j,:)=ddx(j,:)-(w0*(2*i-1))^2*Harm_parameter_a(i,2*j-1)*cos((2*i-1)*w0*Tdata)-(w0*(2*i-1))^2*Harm_parameter_a(i,2*j)*sin((2*i-1)*w0*Tdata);
        end
        x(j,:)=x(j,:)+parameter_a(1,2*j-1);
    end
    x=x';
    xmax_2=x(30001:50001,1);
    xmax(ii).xmax=getmax(xmax_2);
    xmin(ii).xmin=getmin(xmax_2);
    hold on;
    QQ=A(ii)*ones(1,length(xmax(ii).xmax));
    p2=plot(QQ,xmax(ii).xmax,'k.','MarkerSize',5);
    hold on;
end
hold on
%%
load 'bifurcation_unstable_period2_N=50_omega_new=0.471_A_1_-0.01_0.mat'
for ii=1:1:7 % 0.94——1
    parameter_a=every_a(ii).parameter_a;
    A(ii)=every_a(ii).A;
    Tdata1=0:0.01:500;
    % w0=parameter_a(1,1);
    Harm_parameter_a=parameter_a(2:end,:);
    %% 计算方程残差
    x=zeros(N_dof,length(Tdata1));
    dx=zeros(N_dof,length(Tdata1));
    ddx=zeros(N_dof,length(Tdata1));
    % 有X,Y两个自由度
    % for i=1:N_harm   % i=1,3,5
    for j=1:N_dof
        for i=1:N_harm   % i=1,3,5
            x(j,:)=x(j,:)+Harm_parameter_a(i,2*j-1)*cos(i*w_new*Tdata1)+Harm_parameter_a(i,2*j)*sin(i*w_new*Tdata1);
%             dx(j,:)=dx(j,:)-w0*(2*i-1)*Harm_parameter_a(i,2*j-1)*sin((2*i-1)*w0*Tdata)+w0*(2*i-1)*Harm_parameter_a(i,2*j)*cos((2*i-1)*w0*Tdata);
%             ddx(j,:)=ddx(j,:)-(w0*(2*i-1))^2*Harm_parameter_a(i,2*j-1)*cos((2*i-1)*w0*Tdata)-(w0*(2*i-1))^2*Harm_parameter_a(i,2*j)*sin((2*i-1)*w0*Tdata);
        end
        x(j,:)=x(j,:)+parameter_a(1,2*j-1);
    end
    x=x';
    xmax_2=x(30001:50001,1);
    xmax(ii).xmax=getmax(xmax_2);
    xmin(ii).xmin=getmin(xmax_2);
    hold on;
    QQ=A(ii)*ones(1,length(xmax(ii).xmax));
    p3=plot(QQ,xmax(ii).xmax,'b.','MarkerSize',4);
    hold on;
end
hold on

%%
load 'bifurcation_unstable_period2_N=50_omega_new=0.471_A_1_-0.01_0.mat'
for ii=8:1:51 % 0.93——0
    parameter_a=every_a(ii).parameter_a;
    A(ii)=every_a(ii).A;
    Tdata1=0:0.01:500;
    % w0=parameter_a(1,1);
    Harm_parameter_a=parameter_a(2:end,:);
    %% 计算方程残差
    x=zeros(N_dof,length(Tdata1));
    dx=zeros(N_dof,length(Tdata1));
    ddx=zeros(N_dof,length(Tdata1));
    % 有X,Y两个自由度
    % for i=1:N_harm   % i=1,3,5
    for j=1:N_dof
        for i=1:N_harm   % i=1,3,5
            x(j,:)=x(j,:)+Harm_parameter_a(i,2*j-1)*cos(i*w_new*Tdata1)+Harm_parameter_a(i,2*j)*sin(i*w_new*Tdata1);
%             dx(j,:)=dx(j,:)-w0*(2*i-1)*Harm_parameter_a(i,2*j-1)*sin((2*i-1)*w0*Tdata)+w0*(2*i-1)*Harm_parameter_a(i,2*j)*cos((2*i-1)*w0*Tdata);
%             ddx(j,:)=ddx(j,:)-(w0*(2*i-1))^2*Harm_parameter_a(i,2*j-1)*cos((2*i-1)*w0*Tdata)-(w0*(2*i-1))^2*Harm_parameter_a(i,2*j)*sin((2*i-1)*w0*Tdata);
        end
        x(j,:)=x(j,:)+parameter_a(1,2*j-1);
    end
    x=x';
    xmax_2=x(30001:50001,1);
    xmax(ii).xmax=getmax(xmax_2);
    xmin(ii).xmin=getmin(xmax_2);
    hold on;
    QQ=A(ii)*ones(1,length(xmax(ii).xmax));
    p4=plot(QQ,xmax(ii).xmax,'r.','MarkerSize',4);
    hold on;
end
hold on

%%
load 'bifurcation_unstable_period2_N=50_omega_new=0.471_A_1_0.01_2.mat'
for ii=48:1:100 % 1.47——2
    parameter_a=every_a(ii).parameter_a;
    A(ii)=every_a(ii).A;
    Tdata1=0:0.01:500;
    % w0=parameter_a(1,1);
    Harm_parameter_a=parameter_a(2:end,:);
    %% 计算方程残差
    x=zeros(N_dof,length(Tdata1));
    dx=zeros(N_dof,length(Tdata1));
    ddx=zeros(N_dof,length(Tdata1));
    % 有X,Y两个自由度
    % for i=1:N_harm   % i=1,3,5
    for j=1:N_dof
        for i=1:N_harm   % i=1,3,5
            x(j,:)=x(j,:)+Harm_parameter_a(i,2*j-1)*cos(i*w_new*Tdata1)+Harm_parameter_a(i,2*j)*sin(i*w_new*Tdata1);
%             dx(j,:)=dx(j,:)-w0*(2*i-1)*Harm_parameter_a(i,2*j-1)*sin((2*i-1)*w0*Tdata)+w0*(2*i-1)*Harm_parameter_a(i,2*j)*cos((2*i-1)*w0*Tdata);
%             ddx(j,:)=ddx(j,:)-(w0*(2*i-1))^2*Harm_parameter_a(i,2*j-1)*cos((2*i-1)*w0*Tdata)-(w0*(2*i-1))^2*Harm_parameter_a(i,2*j)*sin((2*i-1)*w0*Tdata);
        end
        x(j,:)=x(j,:)+parameter_a(1,2*j-1);
    end
    x=x';
    xmax_2=x(30001:50001,1);
    xmax(ii).xmax=getmax(xmax_2);
    xmin(ii).xmin=getmin(xmax_2);
    hold on;
    QQ=A(ii)*ones(1,length(xmax(ii).xmax));
    p5=plot(QQ,xmax(ii).xmax,'r.','MarkerSize',4);
    hold on;
end
hold off

h1=legend([p1 p2 p3 p4 p5],{'$$stable \; period2$$','$$ unstable \; period1 $$','$$stable \; period2$$','$$stable \; period1$$','$$stable \; period1$$'});
set(h1,'Interpreter','latex','FontSize',15);
% xlabel('$ A $','Interpreter','latex');
% ylabel('$ {\rm Extremes \; of} \; x_1 $','Interpreter','latex')
% h1=legend('$$x_1$$','$$x_2$$');
% set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);









