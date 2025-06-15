clc; clear; close all;

n = 6;          % Robot eklem sayısı

% Her bir eklem uzunluğu (m)
link_uzunlugu = [0.5;
                 0.2;
                 0.2;
                 0.2;
                 0.4;
                 0.2;];

% Eklem eksen yönleri (3x6 matris)
link_yonelimi = [0 0 1;
                 0 1 0;
                 1 0 0;
                 0 0 1;
                 1 0 0;
                 0 0 1;];

% Dönme eksenleri (her eklem için [x y z]) full Y ekseninde dönüyor
h_yonelim = [0 1 0;
             0 1 0;
             0 1 0;
             0 1 0;
             0 1 0;
             0 1 0;];


% Eklem tipleri (1: Revolute, 0: Prismatic)
joint_types = [1;1;1;1;1;1;];

% Robot taban konumu
L_base = [0; 0; 0;];

% dinamik paremetreler
m = [1;1;1;1;1;1;]; % Her bir linkin kütlesi 
kutle_mekezi = [0.5;0.5;0.5;0.5;0.5;0.5;];  % Kütle merkezi linkin ortasında 
I_values = repmat([0.01 0 0 0 0.01 0 0 0 0.01], n, 1); %  atalet tensörü

% Kütle merkezi vektörleri (3 x n)
kutle_merkezi_vektoru = (link_uzunlugu .* kutle_mekezi);
kutle_mekezi = (kutle_merkezi_vektoru .* link_yonelimi)';

% değişkenler
base_pos = [0; 0; 0];
V_b = zeros(6,1);
V_b_dot = zeros(6,1);
theta_dot = zeros(n,1);
dtheta = zeros(n,1);
theta_dot_dot = zeros(n,1);
theta_dot_prev = zeros(n,1);
% zamanlama
t = 0;          
dt = 0.01;     
sim_time = 2;
% cizilecek değişkenler
L = (link_uzunlugu .* link_yonelimi)';
torque_values = zeros(sim_time/dt+1, n);
time_values = zeros(sim_time/dt+1, 1);
points = zeros(n+1, 3);                      % Robot kolu noktaları
link_positions = zeros(n+1, 3, sim_time/dt+1);
uc_islevci_adimlama = zeros(sim_time/dt+1, 3);

% hız ve kuvvetlerimiz
V_t = @(t) [0; 0; 0; 0; cos(2*pi*t); 0];
F_t = @(t) [0; 0; 0; 0; 0; 0];

function skew_symmetric = skew(v)
    skew_symmetric = [0,-v(3), v(2);
                     v(3), 0, -v(1);
                     -v(2),v(1), 0];
end

function R = Rodriguez(dtheta,skew_w)
    R = eye(3) + sin(dtheta) * skew_w + (1 - cos(dtheta)) * skew_w ^ 2;
end

function phi_in = phi_calc(v)
    phi_in = [eye(3), zeros(3,3);
              -skew(v), eye(3)];
end


%simülasyon
step = 1;
while t <= sim_time
    L = link_vektor_guncelle(L, dtheta, joint_types, h_yonelim);
    phi = calculate_phi(L);

    phi_t = zeros(6, 6*n);
    phi_t(:, end-5:end) = phi_calc(L(:,end));

    phi_b = zeros(6*n,6);
    phi_b(1:6,1:6) = phi_calc(L_base);
    phi_tb = phi_t * phi * phi_b;

    H = calculate_H(joint_types, h_yonelim);
    J = phi_t * phi * H;
    V = phi * H * theta_dot;

    a = a_calc(n, V, L, V_b); % ivme
    artik_terim = b_calc(I_values, kutle_mekezi, m, V, n); % artık terim
    M = calculate_M(I_values, kutle_mekezi, m, n); % atalet matrisi
    V_dot = phi * (H * theta_dot_dot + a + phi_b * V_b_dot);

    theta_dot = pinv(J) * (V_t(t) - phi_tb * V_b);
    theta_dot_dot = (theta_dot - theta_dot_prev) / dt;
    dtheta = theta_dot * dt;

    F = phi' * (M * V_dot + artik_terim + phi_t' * F_t(t));
    Torque = H' * F;

    if step <= size(torque_values,1)
        torque_values(step,:) = Torque';
        time_values(step) = t;
    end
    
    points(1,:) = base_pos';
    for i = 1:n
        points(i+1,:) = points(i,:) + L(:,i)';
    end
    link_positions(:,:,step) = points;
    uc_islevci_adimlama(step,:) = points(end,:);

    theta_dot_prev = theta_dot;
    t = t + dt;
    step = step + 1;
end


% tork grafikleri
for i = 1:n
    figure(10+i); clf;
    set(gcf, 'Name',['Ekleme ' num2str(i) ' Tork Grafiği'],'NumberTitle','off');
    plot(time_values, torque_values(:,i), 'b', 'LineWidth', 1.5);
    xlabel('Zaman (s)'); ylabel('Tork (Nm)');
    title(['Ekleme ' num2str(i) ' Tork-Zaman Grafiği']);
    grid on;
end

% robot grafiği
figure(20); clf;
set(gcf, 'Name','Robot Simülasyonu','NumberTitle','off');
view(3);
axis equal;
grid on;
xlabel('X Ekseni (m)'); ylabel('Y Ekseni (m)'); zlabel('Z Ekseni (m)');
title('Robotun Simülasyonu');
xlim([-1 1]); ylim([-1 1]); zlim([-0.2 2]);
hold on;

sabit_uzuv = base_pos' + L_base';
for k = 1:step-1
    cla;
    
    plot3([base_pos(1), sabit_uzuv(1)],...
          [base_pos(2), sabit_uzuv(2)],...
          [base_pos(3), sabit_uzuv(3)],...
          'r*-', 'LineWidth', 2);
    
    link_positions(1,:,k) = sabit_uzuv;
    for i = 2:n+1
        % Bağlantı çizgileri
        plot3([link_positions(i-1,1,k), link_positions(i,1,k)],...
              [link_positions(i-1,2,k), link_positions(i,2,k)],...
              [link_positions(i-1,3,k), link_positions(i,3,k)],...
              'r-', 'LineWidth', 2);
        
        % Eklem noktaları
        scatter3(link_positions(i,1,k),...
                 link_positions(i,2,k),...
                 link_positions(i,3,k),...
                 'MarkerEdgeColor','r');
    end
    drawnow;
end
hold off;



function phi = calculate_phi(L)
    [~, n] = size(L);
    phi = zeros(6*n, 6*n);
    for i = 1:n
        idx = (i-1)*6 + (1:6);
        phi(idx, idx) = eye(6);
        if i > 1
            phi(idx, idx-6) = phi_calc(L(:, i-1));
        end
    end
    for i = 1:n
        for j = i-1:-1:1
            idx_i = (i-1)*6 + (1:6);
            idx_j = (j-1)*6 + (1:6);
            idx_j1 = j*6 + (1:6);
            phi(idx_i, idx_j) = phi(idx_i, idx_j1) * phi(idx_j1, idx_j);
        end
    end
end

function H = calculate_H(joint_types, h_yonelim)
    n = length(joint_types);
    H = zeros(6 * n, n);
    for i = 1:n
        if joint_types(i) == 1
            H(6*(i-1)+1:6*i, i) = [zeros(3,1); h_yonelim(i,:)'];
        else
            H(6*(i-1)+1:6*i, i) = [h_yonelim(i,:)'; zeros(3,1)];
        end
    end
end

function a = a_calc(n,V,L,V_b)
    a = zeros(6 * n, 1);
    g = [0 0 9.81]';
    if V_b == zeros(6,1)
        a(1:6) = [0;0;0;g];
    end
    for i=2:n
        w_k_1 = V((i-2)*6+1:(i-2)*6+3);
        w_k = V((i-1)*6+1:(i-1)*6+3);
        a((i-1)*6+1:(i-1)*6+3) = cross(w_k_1,w_k);
        a((i-1)*6+4:6*i) = cross(w_k_1, cross(w_k_1,L(:,i-1)));
    end
end

function artik_terim = b_calc(I_vals, kutle_mekezi, m, V, n)
    artik_terim = zeros(6 * n, 1);
    for i=1:n
        idx = 6*(i-1)+1:6*(i-1)+3;
        idy = 6*(i-1)+4:6*i;
        w_k = V(idx);
        I = reshape(I_vals(i,:),3,3);
        artik_terim(idx) = cross(w_k, I * w_k);
        artik_terim(idy) = m(i) * cross(w_k, cross(w_k,kutle_mekezi(:,i)));
    end
end

function M = calculate_M(I_vals, kutle_mekezi, m, n)
    M = zeros(6*n,6*n);
    for i=1:n
        idx = (i-1)*6 + (1:6);
        I = reshape(I_vals(i,:),3,3);
        M(idx,idx) = [I, m(i) * skew(kutle_mekezi(:,i));
                      m(i) * skew(kutle_mekezi(:,i)), eye(3)*m(i)];
    end
end

function L = link_vektor_guncelle(L, dtheta, joint_types, h_yonelim)
    for i = 1:length(joint_types)
        if joint_types(i) == 1
            skew_w = skew(h_yonelim(i,:));
            R = Rodriguez(dtheta(i),skew_w);
            L(:, i:end) = R * L(:, i:end);
            if i == length(joint_types)
                L(:, i) = R * L(:, end);
            end
        else
            L(:,i) = L(:,i) + dtheta(i);
        end
    end
end
