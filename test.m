clc
clear
close all
%各种初始化
%信号均值
V1 = 2.2;
V2 = 2.5;
V3 = 2.8;
%信号方差
sigma_1 = 1;
sigma_2 = 1.5;
sigma_3 = 1.2;
%代价函数
C_00 = 0;
C_11 = 0;
C_01 = 1;
C_10 = 1;

% 因为要画出代价函数RB随着P0的变化曲线
P0 = linspace(0,1);
P1 = 1 - P0;

RB =ones(1,100)*10000; 
R1 = zeros(1,100);
R2 = zeros(1,100);
R3 = zeros(1,100);

%终止准则，两个RB函数差距不能超过
Criteion = 10e-3;

%计算代价
C_F = P0 * (C_10 - C_00);
C_D = P1 * (C_01 - C_11);
C   = P1 * C_01 + P0 * C_00;

%门限的初始化，有三个传感器，故有三个、初始门限应该是能够随便给的
%初始化三个传感器的门限值
T1 = zeros(1,100);
T2 = zeros(1,100);
T3 = zeros(1,100);

%%
%生成每一个函数的概率分布函数
z1_H1 = makedist('Normal', V1, sigma_1);
z1_H0 = makedist('Normal', 0 , sigma_1);
z2_H1 = makedist('Normal', V2, sigma_2);
z2_H0 = makedist('Normal', 0 , sigma_2);
z3_H1 = makedist('Normal', V3, sigma_3);
z3_H0 = makedist('Normal', 0 , sigma_3);

%%

for trials = 1:100
     while(1)
        %第一步，传感器根据各自门限计算检测概率
        [P_00_1, P_f1, P_m1, P_d1] =  Get_Probility(z1_H0,z1_H1,T1(trials));
        [P_00_2, P_f2, P_m2, P_d2] =  Get_Probility(z2_H0,z2_H1,T2(trials));
        [P_00_3, P_f3, P_m3, P_d3] =  Get_Probility(z3_H0,z2_H1,T3(trials));


        %根据传感器的检测概率和虚警概率，计算融合中心的P(U|H1) P(U|H0) P(U0|U)
        % 先写一个傻瓜版本
        %这个是在H0的情况下的，8种情况，声明数组
        P_U_H0 = zeros(1, 8);

        P_U_H0(1) = P_00_1 * P_00_2 * P_00_3;
        P_U_H0(2) = P_00_1 * P_00_2 * P_f3;
        P_U_H0(3) = P_00_1 * P_f2 * P_00_3;
        P_U_H0(4) = P_00_1 * P_f2 * P_f3;
        P_U_H0(5) = P_f1 * P_00_2 * P_00_3;
        P_U_H0(6) = P_f1 * P_00_2 * P_f3;
        P_U_H0(7) = P_f1 * P_f2 * P_00_3;
        P_U_H0(8) = P_f1 * P_f2 * P_f2;

        %这个是在H1的情况下的
        P_U_H1 = zeros(1, 8);

        P_U_H1(1) = P_m1 * P_m2 * P_m3;
        P_U_H1(2) = P_m1 * P_m2 * P_d3;
        P_U_H1(3)=  P_m1 * P_d2 * P_m3;
        P_U_H1(4) = P_m1 * P_d2 * P_d3;
        P_U_H1(5) = P_d1 * P_m2 * P_m3;
        P_U_H1(6) = P_d1 * P_m2 * P_d3;
        P_U_H1(7) = P_d1 * P_d2 * P_m3;
        P_U_H1(8) = P_d1 * P_d2 * P_d3;
        %验证是否正确可以sum(P_U_H1) = 1？

        %后面计算传感器的贡献时需要用到
        %写每个传感器的 P(~ui|H0)  P(~ui|H1) 4代表有4项00 01 10 11   3代表的是有三个传感器
        P_Ui_H0 = zeros(1,4,3);
        P_Ui_H1 = zeros(1,4,3);
        P_Ui_H0(:,1,1) = P_00_2 * P_00_3; %代表是 00 
        P_Ui_H0(:,2,1) = P_00_2 * P_f3;;
        P_Ui_H0(:,3,1) = P_f2 * P_00_3;
        P_Ui_H0(:,4,1) = P_f2 * P_f3;

        P_Ui_H0(:,1,2) = P_00_1 *  P_00_3;
        P_Ui_H0(:,2,2) = P_00_1 * P_f3;
        P_Ui_H0(:,3,2) = P_f1   * P_00_3;
        P_Ui_H0(:,4,2) = P_f1   * P_f3;

        P_Ui_H0(:,1,3) = P_00_1 *  P_00_2;
        P_Ui_H0(:,2,3) = P_00_1 *  P_f2;
        P_Ui_H0(:,3,3) = P_f1   *  P_00_2;
        P_Ui_H0(:,4,3) = P_f1   *  P_f2;

        %H1情况下
        P_Ui_H1(:,1,1) =  P_m2 * P_m3;
        P_Ui_H1(:,2,1) =  P_m2 * P_d3;
        P_Ui_H1(:,3,1) =  P_d2 * P_m3;
        P_Ui_H1(:,4,1) =  P_d2 * P_d3;

        P_Ui_H1(:,1,2) = P_m1 * P_m3;
        P_Ui_H1(:,2,2) = P_m1 * P_d3;
        P_Ui_H1(:,3,2) = P_d1 * P_m3;
        P_Ui_H1(:,4,2) = P_d1 * P_d3;

        P_Ui_H1(:,1,3) = P_m1 * P_m2;
        P_Ui_H1(:,2,3) = P_m1 * P_d2;
        P_Ui_H1(:,3,3) = P_d1 * P_m2;
        P_Ui_H1(:,4,3) = P_d1 * P_d2;

        %     C_F * P_000_H0 - C_D * P_000_H1
        %肯定要写成一个数组，不然没法做了。。。。明天再写了2019.4.27 0.10
        %P(U0=1|U)
        P_U0_1_U = zeros(1, 8);
        for i=1:length(P_U0_1_U)
            %C_F(trials)代表的是第trials的代价
            if(C_F(trials) * P_U_H0(i) - C_D(trials) * P_U_H1(i) < 0)
                P_U0_1_U(i) = 1;
            else
                P_U0_1_U(i) = 0;
            end
        end

        %计算代价函数RB,如果要把代价函数设置为
        RB_tmp = 0;
        for i=1:length(P_U0_1_U)
           RB_tmp = RB_tmp +  P_U0_1_U(i) * (C_F(trials) * P_U_H0(i) - C_D(trials) * P_U_H1(i)); 
        end
        RB_tmp = RB_tmp + C(trials);

        %判断是否可以终止的准则
          if(abs(RB(trials) - RB_tmp) > Criteion )            
                %第三步，计算传感器的贡献A（ui）
                A_sensor1 = zeros(1, 4);%传感器的贡献， A(1)表示第一个传感器
                A_sensor2 = zeros(1, 4);
                A_sensor3 = zeros(1, 4);

                A_sensor1(1) = P_U0_1_U(5) - P_U0_1_U(1);
                A_sensor1(2) = P_U0_1_U(6) - P_U0_1_U(2);
                A_sensor1(3) = P_U0_1_U(7) - P_U0_1_U(3);
                A_sensor1(4) = P_U0_1_U(8) - P_U0_1_U(4);


                A_sensor2(1) = P_U0_1_U(3) - P_U0_1_U(1);
                A_sensor2(2) = P_U0_1_U(4) - P_U0_1_U(2);
                A_sensor2(3) = P_U0_1_U(7) - P_U0_1_U(3);
                A_sensor2(4) = P_U0_1_U(8) - P_U0_1_U(4);

                A_sensor3(1) = P_U0_1_U(2) - P_U0_1_U(1);
                A_sensor3(2) = P_U0_1_U(4) - P_U0_1_U(3);
                A_sensor3(3) = P_U0_1_U(6) - P_U0_1_U(5);
                A_sensor3(4) = P_U0_1_U(8) - P_U0_1_U(7);

                A = zeros(1,4,3);
                A(:,:,1) = A_sensor1;
                A(:,:,2) = A_sensor2;
                A(:,:,3) = A_sensor3;

                %计算新的门限
                P_ui_H0 = zeros(1, 4);
                P_ui_H1 = zeros(1, 4);

                Molecule    = zeros(1,3);%每个传感器的分母项
                Denominator = zeros(1,3);
                for i = 1:3 %取决于传感器的个数
                    for j = 1:4 %每一项都等于一个全概率公式
                        Molecule(i)    = Molecule(i)    +  A(:,j,i) * P_Ui_H0(:,j,i);
                        Denominator(i) = Denominator(i) +  A(:,j,i) * P_Ui_H1(:,j,i); 
                    end
                        Molecule(i)    = C_F(trials) * Molecule(i);
                        Denominator(i) = C_D(trials) * Denominator(i);
                end

                T1(trials) = Molecule(1) / Denominator(1);
                T2(trials) = Molecule(2) / Denominator(2);
                T3(trials) = Molecule(3) / Denominator(3);
                RB(trials) = RB_tmp;                                        
          else
              RB(trials) = RB_tmp;
              break;
          end             
     end 
     sprintf("trials = %f",trials)
end
figure(1)
plot(P0,RB,'r'); 
title('贝叶斯风险曲线图');
xlabel('p0');
ylabel('贝叶斯风险值');
legend('融合系统贝叶斯风险曲线');
figure(2)
plot(P0,T1,'r',P0,T2,'g',P0,T3,'b'); 
title('最佳门限');
xlabel('p0');
ylabel('门限值');
legend('T1','T2','T3');