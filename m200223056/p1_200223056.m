clc;
clear;

n1=6; % ilk robot eklem sayosı
n2=6; % 2. robot eklem sayısı

%1. robotun link uzunlukları
initial_L1=[
0.2,0.15,0,0.15,0,0.1;
0,0,0,0,0,0;
0.2,0,0.15,0,0.1,0.1;
];

current_L1=initial_L1;
drawn_L1=transpose(current_L1)

%1. robotun link uzunlukları
initial_L2=[0, -0.15, 0, -0.15, 0, -0.1;
            0, 0,     0, 0,     0, 0;
            0.2, 0, 0.15, 0, 0.1, 0.1];
current_L2=initial_L2;
drawn_L2=transpose(current_L2)

% dönme eksenleri y ekseninde dönüyor hepsi
h1_type = {'rev','rev','rev','rev','rev','rev'}; 
h1 = [0,0,0,0,0,0;
      1,1,0,1,1,1;  
      0,0,1,0,0,0]; 

h2_type = {'rev','rev','rev','rev','rev','rev'}; 
h2 = [0,0,0,0,0,0;
      1,1,0,1,1,1;  
      0,0,1,0,0,0]; 

I=eye(3);
teta1=zeros(n1,1); %teta1 matrisi
teta2=zeros(n2,1); %%teta2 matrisi

% end efektör noktası
ee1_pos = drawn_L1(end,:);
ee2_pos = drawn_L2(end,:);

for t=0:0.01:2
    % önceki adımı silme
    clf;
    hold on;
    
    % Vt linner ve açısal hızlarımız
    V_t1 = [0; 0; 0; 0.9*sin(2*pi*t); 0; 0.9*cos(2*pi*t)];
    V_t2 = [0; 0; 0; 0.9*sin(2*pi*t + pi/4); 0; 0.9*cos(2*pi*t + pi/4)];
    
    % Reset drawing for Robot 1
    current_L1 = initial_L1;
    drawn_L1 = transpose(current_L1);
    
    for i = 1:n1
        if strcmp(h1_type{i}, 'rev')
            teta1(i) = teta1(i) + 0.01 * V_t1(4);
            R = I + sin(teta1(i)) * skewSymmetricMatrix(h1(:,i)) + ...
                    (1 - cos(teta1(i))) * skewSymmetricMatrix(h1(:,i))^2;
            current_L1(:,i) = R * initial_L1(:,i);
        elseif strcmp(h1_type{i}, 'pri')
            teta1(i) = teta1(i) + 0.01 * V_t1(4);
            current_L1(:,i) = initial_L1(:,i) + teta1(i) * h1(:,i);
        end
        drawn_L1(i+1,:) = drawn_L1(i,:) + transpose(current_L1(:,i));
    end

    
    current_L2 = initial_L2;
    drawn_L2 = transpose(current_L2);
    
    
    % rodriges
    for i = 1:n2
        if strcmp(h2_type{i}, 'rev')
            teta2(i) = teta2(i) + 0.01 * V_t2(4);
            R = I + sin(teta2(i)) * skewSymmetricMatrix(h2(:,i)) + ...
                    (1 - cos(teta2(i))) * skewSymmetricMatrix(h2(:,i))^2;
            current_L2(:,i) = R * initial_L2(:,i);
        elseif strcmp(h2_type{i}, 'pri')
            teta1(i) = teta1(i) + 0.01 * V_t2(4);
            current_L2(:,i) = initial_L2(:,i) + teta1(i) * h2(:,i);
        end
        drawn_L2(i+1,:) = drawn_L2(i,:) + transpose(current_L2(:,i));
    end

    
    % end efektörü güncelle
    ee1_pos = drawn_L1(end,:);
    ee2_pos = drawn_L2(end,:);
    
    %cizim
    plot3(drawn_L1(:,1), drawn_L1(:,2), drawn_L1(:,3), 'r-o', 'LineWidth', 2);
    plot3(drawn_L2(:,1), drawn_L2(:,2), drawn_L2(:,3), 'b-o', 'LineWidth', 2);
    
    % iki robot arasında sankiyük taşıyormuş gibi olan bir çizgi
    plot3([ee1_pos(1), ee2_pos(1)], ...
          [ee1_pos(2), ee2_pos(2)], ...
          [ee1_pos(3), ee2_pos(3)], 'k-', 'LineWidth', 3);
    
    % görselleştirme
    view(45,30);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    axis equal;
    xlim([-1 1]); ylim([-1 1]); zlim([0 1.5]);
    grid on;
    title(['sabit platforda çift kollu robot kinematik analizi- Time: ', num2str(t,2)]);
    drawnow;
end

function S = skewSymmetricMatrix(v)
    S = [0, -v(3), v(2); 
         v(3), 0, -v(1); 
         -v(2), v(1), 0];
end