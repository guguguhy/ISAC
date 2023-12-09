close all
clear
clc



%% initialize H
L_n = 3;

N_BS = 128;     % 128 antennas

% 信道增益 cn(0,1)   cn(0,0.01)
g_1 = [GD(0,0.5),GD(0,0.5),GD(0,0.5)];
g_2 = [GD(0,0.05),GD(0,0.05),GD(0,0.05)];
g_3 = [GD(0,0.05),GD(0,0.05),GD(0,0.05)];


% set aoa
theta_1 = SAGL(-70);
theta_2 = SAGL(-40);
theta_3 = SAGL(-10);

A_1 = zeros(N_BS,1);
A_2 = zeros(N_BS,1);
A_3 = zeros(N_BS,1);


for m=1:1:N_BS
    A_1(m) = 1 / sqrt(N_BS) * exp(1j*(m-1)*theta_1);
    A_2(m) = 1 / sqrt(N_BS) * exp(1j*(m-1)*theta_2);
    A_3(m) = 1 / sqrt(N_BS) * exp(1j*(m-1)*theta_3);
end



% BS和第n个用户之间的信道
h_1 = sqrt(N_BS / L_n) * (g_1(1) * A_1' + g_2(1) * A_1' + g_3(1) * A_1');
h_2 = sqrt(N_BS / L_n) * (g_1(2) * A_2' + g_2(2) * A_2' + g_3(2) * A_2');
h_3 = sqrt(N_BS / L_n) * (g_1(3) * A_3' + g_2(3) * A_3' + g_3(3) * A_3');

% 信道矩阵
H = [h_1.',h_2.',h_3.'].';

%% D

M = 400;    % sample points
D = eye(M);


%% Phi

omega = zeros(M,1);
for x=1:1:M
    omega(x) = -1 + (2*x-1)/M;   % sin(x)
end


Phi_trans = zeros(N_BS,M);
for n = 1:1:N_BS
    for m = 1:1:M
        Phi_trans(n,m) = 1 / sqrt(N_BS) * exp(1j*(n-1)*pi*omega(m));
    end
end

Phi = Phi_trans.';

%% b

P_T = 20;
N_RF = 3;       % 3 RF chains
result = sqrt(2*N_RF*P_T/(SAGL(30)-SAGL(10)+SAGL(60)-SAGL(40)));


t_ang = [SAGL(10);SAGL(30);SAGL(40);SAGL(60)];  % target angle

b = zeros(M,1);

% omega_1,omega_2,……,omega_M



for y=1:1:M
    if (omega(y)>= t_ang(1) && omega(y) <= t_ang(2)) || (omega(y)>= t_ang(3) && omega(y) <= t_ang(4))
        b(y) = result;
    else
        b(y) = 0;
    end

end

%% epsilon
epsilon = 0 ;

%% gamma1,2,3(N_c)


f_1 = rand(N_BS,1) + 1j*zeros(N_BS,1);
f_2 = rand(N_BS,1) + 1j*zeros(N_BS,1);
f_3 = rand(N_BS,1) + 1j*zeros(N_BS,1);


gamma_1 = power(abs(h_1 * f_1),2) / (power(abs(h_2 * f_2),2) + power(abs(h_3 * f_3),2)) + epsilon^2;
gamma_2 = power(abs(h_2 * f_2),2) / (power(abs(h_1 * f_1),2) + power(abs(h_3 * f_3),2)) + epsilon^2;
gamma_3 = power(abs(h_3 * f_3),2) / (power(abs(h_1 * f_1),2) + power(abs(h_2 * f_2),2)) + epsilon^2;


%% p

p = exp(1j * linspace(0, 2*pi, M)).';

%% S and S_i

I_nbs = eye(N_BS);
S =  repmat(I_nbs,1,N_RF);


S_s_1 = zeros(N_BS, N_BS*N_RF); % S_1
S_s_2 = zeros(N_BS, N_BS*N_RF); % S_2
S_s_3 = zeros(N_BS, N_BS*N_RF); % S_3

for x=1:1:3
    col_index = (x-1)*N_BS+1 : x*N_BS;
    if x==1
        S_s_1(:,col_index) = eye(N_BS);
    elseif x==2
        S_s_2(:,col_index) = eye(N_BS);
    elseif x==3
        S_s_3(:,col_index) = eye(N_BS);
    end
end

%%
F = [f_1,f_2,f_3];
f = [f_1.',f_2.',f_3.'].';

%%

% A = zeros(M,M);
A = diag(b);

%%



Tou = [30;30;30];


%%
% 通过while循环迭代，直到满足某个条件
epoch = 1;
while epoch <= 1
        
    t = [h_1*S_s_1*f,h_2*S_s_2*f,h_3*S_s_3*f,epsilon].';

    cvx_begin

    num = N_RF*N_BS;

    variable f(num) complex
    minimize(norm(D * (Phi * S * f - A * exp(1j * p)) , 2))

    subject to

    % sum_square(f) <= 20;
    norm(f, 2) <= sqrt(20);
    norm(t,2) <= sqrt(1+1/Tou(1))*t(1);
    norm(t,2) <= sqrt(1+1/Tou(2))*t(2);
    norm(t,2) <= sqrt(1+1/Tou(3))*t(3);

    cvx_end
    epoch = epoch + 1;


    f1_opt = f(1:128);
    f2_opt = f(129:256);
    f3_opt = f(257:384);


    W = pinv(A'*(D'*D)*A)*A'*(D'*D)*Phi;

    p_opt = W * (f1_opt + f2_opt + f3_opt);

    p_fiopt = zeros(400,1);

    for pm=1:1:M
        if abs(p_opt(pm)) ~= 0
            p_fiopt(pm) = p_opt(pm) / abs(p_opt(pm));
        else
            p_fiopt(pm) = 0;
        end
    end

    p = p_fiopt;

end

%%

R = Phi * S * f;

abc = linspace(-90,90,400);
plot(abc,abs(R))


%%

function sin_angle = SAGL(angle)
angle_in_degrees = angle;
angle_in_radians = deg2rad(angle_in_degrees);
sin_angle = sin(angle_in_radians);
end

%%

function complex_samples = GD(mu,variance)
real_part = normrnd(mu, sqrt(variance));
imag_part = normrnd(mu, sqrt(variance));
complex_samples = complex(real_part, imag_part);
end


