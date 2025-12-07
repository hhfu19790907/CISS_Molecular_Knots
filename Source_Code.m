clc;clear;

ti=1;
lambda_SOC=0.36875;
epsilon_n=0;
chiral=1;
R=2;r=1;
p=1;q=2;m=3;
phi_x=0;phi_y=0;phi_z=0;
psi_x=0;psi_y=0;
M=1200;
N_leadL=7;
N_leadR1=15;
N_leadR2=36;
N_leadR3=56;
x0=0;y0=pi/2;z0=pi/2;
c=3e8;


epsilon_0 = 0;
t0 = 1;
t_c = 1;
VL = 1;
VR = 0;
kT = 1;
e = 1.602176634e-19;
h = 4.135667696e-15;
Gamma_L = 1;
Gamma_R = 1;
Gamma_d = 0.005;
df=-1;
num_test=12;
thrshd = 0;






x =@(theta) sin(p*theta+phi_x)+R*sin(q*theta+psi_x);
y =@(theta) cos(p*theta+phi_y)-R*cos(q*theta+psi_y);
z =@(theta) chiral*r*sin(m*theta+phi_z);

theta = linspace(0, 2*pi, M);


dx_dt =@(theta) p*cos(p*theta+phi_x)+q*R*cos(q*theta+psi_x);
dy_dt =@(theta) -p*sin(p*theta+phi_y)+q*R*sin(q*theta+psi_y);
dz_dt =@(theta) chiral*r*m*cos(m*theta+phi_z);


arc_length = trapz(theta, sqrt(dx_dt(theta).^2 + dy_dt(theta).^2+dz_dt(theta) .^2));


N=60;
arc_lengths = linspace(0, arc_length, N+1);
theta_dot=zeros(N+1,1);


for i = 1:N+1

    arc_function =@(theta) trapz(linspace(0, theta, M), sqrt(dx_dt(linspace(0, theta, M)).^2 + dy_dt(linspace(0, theta, M)).^2++ dz_dt(linspace(0, theta, M)).^2));
    theta_dot(i) = fzero(@(theta) arc_function(theta) - arc_lengths(i), [0, 2*pi]);  % 解方程
end

x_print=x(theta);
y_print=y(theta);
z_print=z(theta);


x_dot=x(theta_dot);
y_dot=y(theta_dot);
z_dot=z(theta_dot);

dx_dt1=dx_dt(theta);
dy_dt1=dy_dt(theta);
dz_dt1=dz_dt(theta);



d2x_dt2 =-p^2*sin(p*theta+phi_x)-q^2*R*sin(q*theta+psi_x);
d2y_dt2 =-p^2*cos(p*theta+phi_y)+q^2*R*cos(q*theta+psi_y);
d2z_dt2 =-chiral*r*m^2*sin(m*theta+phi_z);


r_prime=zeros(M,3);
r_double_prime=zeros(M,3);
for i=1:M
    r_prime(i,:) = [dx_dt1(1,i), dy_dt1(1,i), dz_dt1(1,i)];
    r_double_prime(i,:) = [d2x_dt2(i), d2y_dt2(i), d2z_dt2(i)];
end



cross_product = cross(r_prime, r_double_prime);

numerator = sqrt(sum(cross_product.^2, 2));



denominator = sum(r_prime.^2, 2).^(3/2);


epsilon = 1e-10;
kappa = numerator ./ (denominator + epsilon);

kappa_1=zeros(N,1);
for i=1:N
    kappa_1(i)=kappa(round(theta_dot(i)*M/(2*pi))+1);
end

T=r_prime;
T_mag = sqrt(sum(T.^2, 2));
T_unit = T ./ T_mag;


N1=r_double_prime;
N1_mag = sqrt(sum(N.^2, 2));
N1_unit = N1 ./ N1_mag;


B = cross(T_unit, N1_unit);




sigma_1=[1 0; 0 -1];
sigma_2=[0 -1i; 1i, 0];
sigma_3=[0 -1; -1 0];




B_1=zeros(N,3);
sigma_B=cell(1,N);
for i=1:N
    B_1(i,1)=B(round(theta_dot(i)*M/(2*pi))+1,1);
    B_1(i,2)=B(round(theta_dot(i)*M/(2*pi))+1,2);
    B_1(i,3)=B(round(theta_dot(i)*M/(2*pi))+1,3);
    sigma_B{1,i}=B_1(i,1)*sigma_1+B_1(i,2)*sigma_2+B_1(i,3)*sigma_3;
end


Hamiltonian_1=zeros(2*N,2*N);
sigma_D=zeros(2,2);
for i=1:N-1
    sigma_D=(sigma_B{1,i}+sigma_B{1,i+1})*kappa_1(i+1,1);
    Hamiltonian_1(2*i-1,2*i+1)=sigma_D(1,1);
    Hamiltonian_1(2*i-1,2*i+2)=sigma_D(1,2);
    Hamiltonian_1(2*i,2*i+1)=sigma_D(2,1);
    Hamiltonian_1(2*i,2*i+2)=sigma_D(2,2);
end

sigma_D=(sigma_B{1,1}+sigma_B{1,N})*kappa_1(1,1);
Hamiltonian_1(1,2*N-1)=sigma_D(1,1);
Hamiltonian_1(1,2*N)=sigma_D(1,2);
Hamiltonian_1(2,2*N-1)=sigma_D(2,1);
Hamiltonian_1(2,2*N)=sigma_D(2,2);


for i=1:N-1
    sigma_D=-(sigma_B{1,i}+sigma_B{1,i+1})*kappa_1(i,1);
    Hamiltonian_1(2*i+1,2*i-1)=sigma_D(1,1);
    Hamiltonian_1(2*i+1,2*i)=sigma_D(1,2);
    Hamiltonian_1(2*i+2,2*i-1)=sigma_D(2,1);
    Hamiltonian_1(2*i+2,2*i)=sigma_D(2,2);
end

sigma_D=-(sigma_B{1,1}+sigma_B{1,N})*kappa_1(N,1);
Hamiltonian_1(2*N-1,1)=sigma_D(1,1);
Hamiltonian_1(2*N,1)=sigma_D(1,2);
Hamiltonian_1(2*N-1,2)=sigma_D(2,1);
Hamiltonian_1(2*N,2)=sigma_D(2,2);






Hamiltonian_t=zeros(2*N,2*N);
for i=1:N-1
    Hamiltonian_t(2*i-1,2*i+1)=ti;
    Hamiltonian_t(2*i,2*i+2)=ti;
end

for i=1:N-1
    Hamiltonian_t(2*i+1,2*i-1)=ti;
    Hamiltonian_t(2*i+2,2*i)=ti;
end

Hamiltonian_t(1,2*N-1)=ti;
Hamiltonian_t(2,2*N)=ti;
Hamiltonian_t(2*N-1,1)=ti;
Hamiltonian_t(2*N,2)=ti;





Hamiltonian_SOC_half=lambda_SOC*1i*Hamiltonian_1;

Hamiltonian_SOC=Hamiltonian_SOC_half+Hamiltonian_SOC_half';

Hamiltonian_C=zeros(2*N,2*N);
for i=1:2*N
    Hamiltonian_C(i,i)=epsilon_n;
end

Hamiltonian=Hamiltonian_C+Hamiltonian_t+Hamiltonian_SOC;
[V,D]=eig(Hamiltonian);
m=zeros(2*N,1);
for i=1:2*N
    m(i,1)=i;
end
E_D = diag(D);


E_min = D(1,1); E_max = D(2*N,2*N);

E = linspace(E_min,E_max,1000);


Gamma_D_diag = sparse(Gamma_d*eye(2*N));


Gamma_D = zeros(2*N,2*N,2*N);
for ii = 1:2*N
    Gamma_D(ii, ii, ii) = Gamma_d;
end


Gamma_self_LR = zeros(2*N,2*N);
Gamma_self_LR(2*N_leadL-1,2*N_leadL-1) = Gamma_L;
Gamma_self_LR(2*N_leadL,2*N_leadL) = Gamma_L;
Gamma_self_LR(2*N_leadR1-1,2*N_leadR1-1) = Gamma_R;
Gamma_self_LR(2*N_leadR1,2*N_leadR1) = Gamma_R;
Gamma_self_LR(2*N_leadR2-1,2*N_leadR2-1) = Gamma_R;
Gamma_self_LR(2*N_leadR2,2*N_leadR2) = Gamma_R;
Gamma_self_LR(2*N_leadR3-1,2*N_leadR3-1) = Gamma_R;
Gamma_self_LR(2*N_leadR3,2*N_leadR3) = Gamma_R;



Sigma = -1i/2*(Gamma_D_diag + Gamma_self_LR);


G_L_up=zeros(2*N,2*N);
G_L_dn=zeros(2*N,2*N);
G_L_up(2*N_leadL-1,2*N_leadL-1)=Gamma_L;
G_L_dn(2*N_leadL,2*N_leadL)=Gamma_L;
G_L_up=sparse(G_L_up);
G_L_dn=sparse(G_L_dn);
G_L = G_L_up + G_L_dn;

G_R_up=zeros(2*N,2*N);
G_R_down=zeros(2*N,2*N);
G_R_up((2*N_leadR1-1),(2*N_leadR1-1))=Gamma_R;
G_R_up((2*N_leadR2-1),(2*N_leadR2-1))=Gamma_R;
G_R_up((2*N_leadR3-1),(2*N_leadR3-1))=Gamma_R;
G_R_down((2*N_leadR1),(2*N_leadR1))=Gamma_R;
G_R_down((2*N_leadR2),(2*N_leadR2))=Gamma_R;
G_R_down((2*N_leadR3),(2*N_leadR3))=Gamma_R;
G_R_up=sparse(G_R_up);
G_R_down=sparse(G_R_down);
G_R = G_R_up + G_R_down;
G_LR = G_L + G_R;


Con_up = zeros(1, length(E));
Con_down = zeros(1, length(E));
P = zeros(1, length(E));
T_e = zeros(2*N,2*N);
I_test = zeros(1,length(E));
num_test = 3;


identity_2N = speye(2*N);

for e = 1:length(E)

    z1 = E(e) + 1i * 1e-10;
    G = (z1 * identity_2N - Hamiltonian - Sigma) \ identity_2N;
    Ga = G';


    T_e = zeros(2*N,2*N);
    for ii = 1:2*N
        for jj = 1:2*N
            T_e(ii,jj) = Gamma_d^2*abs(G(ii,jj))^2;
        end
    end


    if Gamma_d == 0
        Vn = zeros(N,1);
    else
        M1 = zeros(N,N);
        for ii = 1:N
            for jj = 1:N
                if ii == jj
                    M1(ii,jj) = sum(sum(T_e([2*ii-1,2*ii],setdiff(1:2*N,[2*ii-1,2*ii]))))...
                        + trace(sparse(squeeze(Gamma_D(2*ii-1,:,:)+Gamma_D(2*ii,:,:)))*G*G_LR*Ga);
                else
                    M1(ii,jj) = -sum(sum(T_e([2*ii-1,2*ii],[2*jj-1,2*jj])));
                end
            end
        end
        B1 = zeros(N,1);
        for ii = 1:N
            B1(ii) = trace(sparse(squeeze(Gamma_D(2*ii-1,:,:)+Gamma_D(2*ii,:,:)))*G*G_L*Ga)*VL...
                + trace(sparse(squeeze(Gamma_D(2*ii-1,:,:)+Gamma_D(2*ii,:,:)))*G*G_R*Ga)*VR;
        end
        Vn = M1\B1;
    end


    Con_up(e) = 0;
    for ii = 1:N
        Con_up(e) = Con_up(e) + trace(G_R_up*G*sparse(squeeze(Gamma_D(2*ii-1,:,:)+Gamma_D(2*ii,:,:)))*Ga)*Vn(ii);
    end
    Con_up(e) = Con_up(e) + trace(G_R_up*G*G_L*Ga)*VL;
    Con_up(e) = real(Con_up(e)/(VL-VR));

    Con_down(e) = 0;
    for ii = 1:N
        Con_down(e) = Con_down(e) + trace(G_R_down*G*sparse(squeeze(Gamma_D(2*ii-1,:,:)+Gamma_D(2*ii,:,:)))*Ga)*Vn(ii);
    end
    Con_down(e) = Con_down(e) + trace(G_R_down*G*G_L*Ga)*VL;
    Con_down(e) = real(Con_down(e)/(VL-VR));

    if abs(Con_up(e)) < thrshd
        Con_up(e) = 0;
    end
    if abs(Con_down(e)) < thrshd
        Con_down(e) = 0;
    end

    if (Con_up(e) + Con_down(e)) ~= 0
        P(e) = real((Con_up(e) - Con_down(e)) / (Con_up(e) + Con_down(e)));
    else
        P(e) = 0;
    end


end







% 可视化----------------------------------------------------------------
%% 结构示意图
figure;
hold on;
plot3(x_dot, y_dot, z_dot, 'r.','MarkerFaceColor','r','MarkerSize',15);
plot3(x_print, y_print, z_print, 'b', 'LineWidth', 2);        % 绘制三叶结
hold off;

xlabel('X'); ylabel('Y'); zlabel('Z');
title('Trefoil Knot with Right-Handed Helix');
grid on;
view(3);
axis equal;
rotate3d on;


figure;
quiver3(x_print, y_print, z_print, T_unit(:,1).', T_unit(:,2).', T_unit(:,3).', 0.5, 'r'); hold on;
quiver3(x_print, y_print, z_print, N1_unit(:,1).', N1_unit(:,2).', N1_unit(:,3).', 0.5, 'g');
quiver3(x_print, y_print, z_print, B(:,1).', B(:,2).', B(:,3).', 0.5, 'b');
plot3(x_print, y_print, z_print, 'k'); axis equal;
xlabel('x'); ylabel('y'); zlabel('z');
title('Frenet-Serret Frame Along Variable-Pitch Helix');
legend('Tangent', 'Normal', 'Binormal');



% %% 能级图1
% clc;
% figure;
% plot(zeros(2*N,1),E_D,'b.');
% xlabel(' '); ylabel('energy');
% title('能级分布情况（合）');
%
% %% 能级图2
% clc;
% figure;
% plot(m,E_D,'b.');
% xlabel(' '); ylabel('energy');
% title('能级分布情况（分立）');


%% 电导和极化率
figure;

h1 = plot(E, Con_up, 'k-', 'LineWidth', 1.5); % 上自旋电导
hold on;
h2 = plot(E, Con_down, 'b-', 'LineWidth', 1.5); % 下自旋电导
ylabel('约化电导或自旋极化率'); % 左 y 轴标签

figure;
h3 = plot(E, P, 'r-', 'LineWidth', 1.5); % 自旋极化度

% 添加标题和 x 轴标签
title('自旋极化与能量的关系');
xlabel('能量');

% 显示网格线
grid on;

% 添加图例
legend([h1, h2, h3], {'Con_{up}', 'Con_{dn}', 'P'});

%% 画出曲率随角度的变化
figure;
plot(theta_dot(1:N,1),kappa_1(:,1),'b-', 'LineWidth', 1.5);
title('三叶结中的曲率变化');
xlabel('\theta');
ylabel('\kappa');
hold on;

