clear all;clc;
tic
H = [1 0 0 0;0 0 1 0];
q=2;Len=1000;
N = Len;T = 1;
MC=200;
p=4;
ex1=zeros(MC,N);ex2=zeros(MC,N);
ey1=zeros(MC,N);ey2=zeros(MC,N);% 滤波误差初始化
%蒙特卡罗仿真
for mn=1:MC
    %% 初始x数据产生
    x0 = [1000,10,1000,-10]';
    p0=[ 100 0 0 0;0 1 0 0;0 0 100 0;0 0 0 1;];
    xA = [];
    zA = [x0(1);x0(3)];
    %model-1,匀速运动
    A1 = [1,T,0,0;
        0,1,0,0;
        0,0,1,T;
        0,0,0,1];
    G1=[T^2/2,    0;
        T,      0;
        0,      T^2/2;
        0,      T] ;
    Q1=[0.1^2 0;
        0 0.1^2];
    %model-2,匀速转弯模型
    omgia1=pi/270;
  %  omgia1=pi/90;
    A2=[1 sin(omgia1)/omgia1 0 (cos(omgia1)-1)/omgia1;
        0 cos(omgia1) 0 -sin(omgia1);
        0 -(cos(omgia1)-1)/omgia1 1 sin(omgia1)/omgia1;
        0 sin(omgia1) 0 cos(omgia1);];
    G2=[1/2 0;
        1 0;
        0 1/2;
        0 1;];
    Q2=[0.0144^2 0;
        0 0.0144^2];
    % 产生真实数据
    x = x0;
    for k = 1:350%匀速直线
        x = A1*x + G1*sqrt(Q1)*[randn,randn]';
        xA =[xA x];
    end
    for k = 1:220%匀速圆周转弯
        x = A2*x + G2*sqrt(Q2)*[randn,randn]';
        xA =[xA x];
    end
    for k = 1:430%匀速直线
        x = A1*x + G1*sqrt(Q1)*[randn,randn]';
        xA =[xA x];
    end
    
    %% 量测噪声产生
%     v1 = randn(q,Len) * 2;
%     v2 = randn(q,Len) * 200;
    v1 = randn(q,Len) * 0.01;
    v2 = randn(q,Len) * 200;
    for ii=1:q
        numrand = rand(1,Len);
        v(ii,:) = (numrand>=0.9).*v2(ii,:) + (numrand<0.9).*v1(ii,:);
    end % 冲击噪声   
%     v = randn(q,Len) * 20;% 高斯噪声
     R=v*v'/Len;
     
  %% 量测产生
    for i=1:Len
       zA(:,i) = H*xA(:,i) +v(:,i);
    end
  %% 初始状态给定 
      %x0_hj=[zA(1,2),zA(1,2)-zA(1,1),zA(2,2),zA(2,2)-zA(2,1)]';
      x0_hj=[500,10,500,-10]';
    %% CV模型卡尔曼滤波
    x0_kf_cv=x0_hj;
    P0_kf_cv=p0;%初始协方差
    for i=1:Len
        [P_kv,P_k_k_1v,X1_kv]=kalman(A1,G1,H,Q1,R,zA(:,i),x0_kf_cv,P0_kf_cv);  %%
        X_kv(:,i)=X1_kv;
        x0_kf_cv=X1_kv;P0_kf_cv=P_kv;
    end
    %% MCKFCV 单模型滤波
    x0_mckf_cv = x0_hj;
    P0_mckf_cv = p0;
    for i=1:Len
        [Pkk_MCKF3_CV,X_MCKF3_CV]=MCKF3(A1,G1,H,Q1,R,zA(:,i),x0_mckf_cv,P0_mckf_cv,p);
        X_MCKF_CV(:,i)=X_MCKF3_CV;
        P0_mckf_cv=Pkk_MCKF3_CV;
        x0_mckf_cv=X_MCKF3_CV;
    end
    
    %% CT模型卡尔曼滤波
    omgia1=pi/180;
  %   omgia1=pi/60;
    A2=[1 sin(omgia1)/omgia1 0 (cos(omgia1)-1)/omgia1;
        0 cos(omgia1) 0 -sin(omgia1);
        0 -(cos(omgia1)-1)/omgia1 1 sin(omgia1)/omgia1;
        0 sin(omgia1) 0 cos(omgia1);] ;    
    x0_kf_ct = x0_hj;
    P0_kf_ct = p0; % CT 模型滤波参数
    for i=1:Len
        [P_kt,P_k_k_1t,X1_kt]=kalman(A2,G2,H,Q2,R,zA(:,i),x0_kf_ct,P0_kf_ct);
        X_kt(:,i)=X1_kt;
        x0_kf_ct=X1_kt;P0_kf_ct=P_kt;
    end
    
    %% MCKFCT 单模型滤波
    x0_mckf_ct = x0_hj;
    P0_mckf_ct = p0;
    for i=1:Len
        [Pkk_MCKF3_CT,X_MCKF3_CT]=MCKF3(A2,G2,H,Q2,R,zA(:,i),x0_mckf_ct,P0_mckf_ct,p);
        X_MCKF_CT(:,i)=X_MCKF3_CT;
        P0_mckf_ct=Pkk_MCKF3_CT;
        x0_mckf_ct=X_MCKF3_CT;
    end
    
    ex_kf_cv(mn,:)=X_kv(1,:)-xA(1,:);ey_kf_cv(mn,:)=X_kv(3,:)-xA(3,:);
    ex_kf_ct(mn,:)=X_kt(1,:)-xA(1,:);ey_kf_ct(mn,:)=X_kt(3,:)-xA(3,:);
    ex_mckf_cv(mn,:)=X_MCKF_CV(1,:)-xA(1,:);ey_mckf_cv(mn,:)=X_MCKF_CV(3,:)-xA(3,:);
    ex_mckf_ct(mn,:)=X_MCKF_CT(1,:)-xA(1,:);ey_mckf_ct(mn,:)=X_MCKF_CT(3,:)-xA(3,:);
    
   %% 转移矩阵
    Pi =[0.99,0.01;0.01,0.99];%转移概率
   %  Pi =[0.95,0.05;0.05,0.95];
   %  Pi =[0.9,0.1;0.1,0.9];
      
    %% %%%%%%%%IMM滤波初始化参数%%%%%%% %%%%
    u1=1/2;
    u2=1/2;%2个模型间初始概率传播参数
    x0_imm =x0_hj;
    x1_k_1=x0_imm;x2_k_1=x0_imm; %模型的状态传播参数
    P0_imm = p0;
    P1=P0_imm;P2=P0_imm; %每个模型的状态传播参数
    for k = 1:Len%1-400秒
        %混合概率
        c1=Pi(1,1)*u1+Pi(2,1)*u2;
        c2=Pi(1,2)*u1+Pi(2,2)*u2;
        u11=Pi(1,1)*u1/c1;u12=Pi(1,2)*u1/c2;
        u21=Pi(2,1)*u2/c1;u22=Pi(2,2)*u2/c2;
        x1_m = x1_k_1*u11+x2_k_1*u21;
        x2_m = x1_k_1*u12+x2_k_1*u22;
        p1_k_1=(P1+(x1_k_1-x1_m)*(x1_k_1-x1_m)')*u11+(P2+(x2_k_1-x1_m)*(x2_k_1-x1_m)')*u21;
        p2_k_1=(P1+(x1_k_1-x2_m)*(x1_k_1-x2_m)')*u12+(P2+(x2_k_1-x2_m)*(x2_k_1-x2_m)')*u22;
        %状态预测
        x1_pk1=A1*x1_m; x2_pk1=A2*x2_m;
        p1_k=A1*p1_k_1*A1'+G1*Q1*G1';
        p2_k=A2*p2_k_1*A2'+G2*Q2*G2';
        %预测残差及协方差计算
        zk=zA(:,k);
        v1=zk-H*x1_pk1; v2=zk-H*x2_pk1;
        Sv1=H*p1_k*H'+R;Sv2=H*p2_k*H'+R;
        like1=det(2*pi*Sv1)^(-0.5)*exp(-0.5*v1'*inv(Sv1)*v1);
        like2=det(2*pi*Sv2)^(-0.5)*exp(-0.5*v2'*inv(Sv2)*v2);
        %滤波更新
        K1=p1_k*H'*inv(Sv1); K2=p2_k*H'*inv(Sv2);
        xk1=x1_pk1+K1*v1;xk2=x2_pk1+K2*v2;
        P1=p1_k-K1*Sv1*K1';P2=p2_k-K2*Sv2*K2';
        %模型概率更新
        C=like1*c1+like2*c2;
%         temp_u1=u1;temp_u2=u2; %1
        u1=like1*c1/C;u2=like2*c2/C;
        
%         delta_u=u1-temp_u1;%2
%             Pi(2,1)=Pi(2,1)+delta_u*0.2;
%             Pi(2,2)=1-Pi(2,1);
%             Pi(1,2)=Pi(2,1);
%             Pi(1,1)=Pi(2,2);  %3 
%        if(Pi(1,1)>0.99)
%            Pi(1,1)=0.99;
%            Pi(2,2)=Pi(1,1);
%            Pi(2,1)=1-Pi(1,1);
%            Pi(1,2)=Pi(2,1);
%        end    
%         if(Pi(1,1)<0.95)
%            Pi(1,1)=0.95;
%            Pi(2,2)=Pi(1,1);
%            Pi(2,1)=1-Pi(1,1);
%            Pi(1,2)=Pi(2,1);
%        end
        %估计融合
        xk=xk1*u1+xk2*u2;
        %迭代
        x1_k_1=xk1;x2_k_1=xk2;
        X_imm(:,k)=xk;
        um1(k)=u1;um2(k)=u2;
    end
    ex_imm(mn,:)=X_imm(1,:)-xA(1,:);
    ey_imm(mn,:)=X_imm(3,:)-xA(3,:);
    UM1(mn,:)=um1;UM2(mn,:)=um2;
    
    %% %%%%%%%%IMM_MCKF滤波初始化参数%%%%%%% %%%%
    u1_MCKF=1/2;   u2_MCKF=1/2;
    x0_mckf_imm = x0_hj;
    x1_k_1_MCKF=x0_mckf_imm;x2_k_1_MCKF=x0_mckf_imm; %每个模型的状态传播参数
    P0_mckf_imm =p0;
    P1_MCKF=P0_mckf_imm;P2_MCKF=P0_mckf_imm; %每个模型的状态传播参数

    for k = 1:Len%1-400秒
        %混合概率
        c1_MCKF=Pi(1,1)*u1_MCKF+Pi(2,1)*u2_MCKF;
        c2_MCKF=Pi(1,2)*u1_MCKF+Pi(2,2)*u2_MCKF;
        u11_MCKF=Pi(1,1)*u1_MCKF/c1_MCKF;u12_MCKF=Pi(1,2)*u1_MCKF/c2_MCKF;
        u21_MCKF=Pi(2,1)*u2_MCKF/c1_MCKF;u22_MCKF=Pi(2,2)*u2_MCKF/c2_MCKF;
        x1_m_MCKF = x1_k_1_MCKF*u11_MCKF+x2_k_1_MCKF*u21_MCKF;
        x2_m_MCKF = x1_k_1_MCKF*u12_MCKF+x2_k_1_MCKF*u22_MCKF;
        p1_k_1_MCKF=(P1_MCKF+(x1_k_1_MCKF-x1_m_MCKF)*(x1_k_1_MCKF-x1_m_MCKF)')*u11_MCKF+(P2_MCKF+(x2_k_1_MCKF-x1_m_MCKF)*(x2_k_1_MCKF-x1_m_MCKF)')*u21_MCKF;
        p2_k_1_MCKF=(P1_MCKF+(x1_k_1_MCKF-x2_m_MCKF)*(x1_k_1_MCKF-x2_m_MCKF)')*u12_MCKF+(P2_MCKF+(x2_k_1_MCKF-x2_m_MCKF)*(x2_k_1_MCKF-x2_m_MCKF)')*u22_MCKF;
        %状态预测
        x1_pk1_MCKF=A1*x1_m_MCKF; x2_pk1_MCKF=A2*x2_m_MCKF;
        p1_k_MCKF=A1*p1_k_1_MCKF*A1'+G1*Q1*G1';
        p2_k_MCKF=A2*p2_k_1_MCKF*A2'+G2*Q2*G2';
        % 对应于各模型进行DMCKF
        t=1;
        sigma=3;
        epsilo=0.01;
        zk=zA(:,k);
        
        Bpk_1 = chol(p1_k_MCKF)';
        Brk_1 = chol(R)';
        Bkk_1 = blkdiag(Bpk_1,Brk_1);
        dik_1 = inv(Bkk_1)*[x1_pk1_MCKF;zk];
        wik_1 = [inv(Bpk_1);inv(Brk_1)* H];
        
        xkk_1 = x1_pk1_MCKF;
        while(t<1000)
            temp=xkk_1;
            t=t+1;
            ee_1 = dik_1 - wik_1 * xkk_1;
            Gee_1 = exp(-ee_1.^2/2/sigma^2);
           % Gee_1 =(dik_1.*(wik_1 * xkk_1)+10).^2; 
             %Gee_1 =dik_1.*(wik_1 * xkk_1);
           %  Gee_1 =tanh(0.2*dik_1.*(wik_1 * xkk_1)+5);
            % Gee_1 = 0.9*exp(-ee_1.^2/2/sigma^2)+0.1*tanh(0.2*dik_1.*(wik_1 * xkk_1)+5);
             %Gee_1=0.7*tanh(0.2*dik_1.*(wik_1 * xkk_1)+5)+0.3*dik_1.*(wik_1 * xkk_1);
            Cx_1 = diag(Gee_1(1:4));
            Cy_1 = diag(Gee_1(4+1:6));
%             Pke_hat_1 = inv(Bpk_1)' * (Cx_1) * inv(Bpk_1);
%             R_hat_1 = inv(Brk_1)' * (Cy_1) * inv(Brk_1);
%             Gk_MCC_1 = inv(H'*R_hat_1*H+Pke_hat_1)*H'*R_hat_1;
            
        Pke_hat_1 = Bpk_1 * inv(Cx_1) * Bpk_1';
        R_hat_1 = Brk_1 * inv(Cy_1) * Brk_1';
        Gk_MCC_1 = Pke_hat_1 * H' * inv(H * Pke_hat_1 * H'+R_hat_1);
            
            xkk_1 = x1_pk1_MCKF + Gk_MCC_1*(zk-H *x1_pk1_MCKF );
            
            xe_MCKF_1(:,k) = xkk_1;
            
            if(abs(xe_MCKF_1(:,k)-temp)/abs(temp)<=epsilo)
                break;
            else
                t=t-1;
            end
        end
        %% 2
        Bpk_2 = chol(p2_k_MCKF)';
        Brk_2 = chol(R)';
        Bkk_2 = blkdiag(Bpk_2,Brk_2);
        dik_2 = inv(Bkk_2)*[x2_pk1_MCKF;zk];
        wik_2 = [inv(Bpk_2);inv(Brk_2)* H];
        
        xkk_2 = x2_pk1_MCKF;
        while(t<1000)
            temp=xkk_2;
            t=t+1;
            ee_2 = dik_2 - wik_2 * xkk_2;
            Gee_2 = exp(-ee_2.^2/2/sigma^2)
           %Gee_2 =(dik_2.*(wik_2 * xkk_2)+10).^2
            %Gee_2 =dik_2.*(wik_2 * xkk_2)
            %Gee_2 =tanh(0.2*dik_2.*(wik_2 * xkk_2)+5)
          % Gee_2 = 0.9*exp(-ee_2.^2/2/sigma^2)+0.1*tanh(0.2*dik_2.*(wik_2 * xkk_2)+5)
            %Gee_2=0.7*tanh(0.2*dik_2.*(wik_2 * xkk_2)+5)+0.3*dik_2.*(wik_2 * xkk_2)
            Cx_2 = diag(Gee_2(1:4));
            Cy_2 = diag(Gee_2(4+1:6));
%             Pke_hat_2 = inv(Bpk_2)' * (Cx_2) * inv(Bpk_2);
%             R_hat_2 = inv(Brk_2)' * (Cy_2) *inv(Brk_2);
%             Gk_MCC_2 =inv(H'*R_hat_2*H+Pke_hat_2)*H'*R_hat_2;
            
        Pke_hat_2 = Bpk_2 * inv(Cx_2) * Bpk_2';
        R_hat_2 = Brk_2 * inv(Cy_2) * Brk_2';
        Gk_MCC_2 = Pke_hat_2 * H' * inv(H * Pke_hat_2 * H'+R_hat_2);
            
        
            xkk_2 = x2_pk1_MCKF + Gk_MCC_2*(zk-H *x2_pk1_MCKF );
            
            xe_MCKF_2(:,k) = xkk_2;
            
            if(abs(xe_MCKF_2(:,k)-temp)/abs(temp)<=epsilo)
                break;
            else
                t=t-1;
            end
        end
     
        %预测残差及协方差计算
        v1_MCKF=zk-H*x1_pk1_MCKF;
        v2_MCKF=zk-H*x2_pk1_MCKF;
        Sv1_MCKF=H*p1_k_MCKF*H'+R;
        Sv2_MCKF=H*p2_k_MCKF*H'+R;
        like1_MCKF=det(2*pi*Sv1_MCKF)^(-0.5)*exp(-0.5*v1_MCKF'*inv(Sv1_MCKF)*v1_MCKF);
        like2_MCKF=det(2*pi*Sv2_MCKF)^(-0.5)*exp(-0.5*v2_MCKF'*inv(Sv2_MCKF)*v2_MCKF);
        %滤波更新
        K1_MCKF=Gk_MCC_1;K2_MCKF=Gk_MCC_2;
        xk1_MCKF=xe_MCKF_1(:,k);xk2_MCKF=xe_MCKF_2(:,k);
        P1_MCKF=(eye(4) - Gk_MCC_1 * H) * p1_k_MCKF * (eye(4) - Gk_MCC_1 * H)'+ Gk_MCC_1 * R * Gk_MCC_1';
        P2_MCKF=(eye(4) - Gk_MCC_2 * H) * p2_k_MCKF * (eye(4) - Gk_MCC_2 * H)'+ Gk_MCC_2 * R * Gk_MCC_2';
        %模型概率更新
        C_MCKF=like1_MCKF*c1_MCKF+like2_MCKF*c2_MCKF;
        u1_MCKF=like1_MCKF*c1_MCKF/C_MCKF;
        u2_MCKF=like2_MCKF*c2_MCKF/C_MCKF;
        %估计融合
        xk_MCKF=xk1_MCKF*u1_MCKF+xk2_MCKF*u2_MCKF;
        %迭代
        x1_k_1_MCKF=xk1_MCKF;x2_k_1_MCKF=xk2_MCKF;
        X_imm_MCKF(:,k)=xk_MCKF;
        um1_MCKF(k)=u1_MCKF;um2_MCKF(k)=u2_MCKF;
    end
    
      %% %%%%%%%%IMM_SMCKF滤波初始化参数%%%%%%% %%%%
    u1_SMCKF=1/2;   u2_SMCKF=1/2;
    x0_smckf_imm = x0_hj;
    x1_k_1_SMCKF=x0_smckf_imm;x2_k_1_SMCKF=x0_smckf_imm; %每个模型的状态传播参数
    P0_smckf_imm =p0;
    P1_SMCKF=P0_smckf_imm;P2_SMCKF=P0_smckf_imm; %每个模型的状态传播参数

    for k = 1:Len%1-400秒
        %混合概率
        c1_SMCKF=Pi(1,1)*u1_SMCKF+Pi(2,1)*u2_SMCKF;
        c2_SMCKF=Pi(1,2)*u1_SMCKF+Pi(2,2)*u2_SMCKF;
        u11_SMCKF=Pi(1,1)*u1_SMCKF/c1_SMCKF;u12_SMCKF=Pi(1,2)*u1_SMCKF/c2_SMCKF;
        u21_SMCKF=Pi(2,1)*u2_SMCKF/c1_SMCKF;u22_SMCKF=Pi(2,2)*u2_SMCKF/c2_SMCKF;
        x1_m_SMCKF = x1_k_1_SMCKF*u11_SMCKF+x2_k_1_SMCKF*u21_SMCKF;
        x2_m_SMCKF = x1_k_1_SMCKF*u12_SMCKF+x2_k_1_SMCKF*u22_SMCKF;
        p1_k_1_SMCKF=(P1_SMCKF+(x1_k_1_SMCKF-x1_m_SMCKF)*(x1_k_1_SMCKF-x1_m_SMCKF)')*u11_SMCKF+(P2_SMCKF+(x2_k_1_SMCKF-x1_m_SMCKF)*(x2_k_1_SMCKF-x1_m_SMCKF)')*u21_SMCKF;
        p2_k_1_SMCKF=(P1_SMCKF+(x1_k_1_SMCKF-x2_m_SMCKF)*(x1_k_1_SMCKF-x2_m_SMCKF)')*u12_SMCKF+(P2_SMCKF+(x2_k_1_SMCKF-x2_m_SMCKF)*(x2_k_1_SMCKF-x2_m_SMCKF)')*u22_SMCKF;
        %状态预测
        x1_pk1_SMCKF=A1*x1_m_SMCKF; x2_pk1_SMCKF=A2*x2_m_SMCKF;
        p1_k_SMCKF=A1*p1_k_1_SMCKF*A1'+G1*Q1*G1';
        p2_k_SMCKF=A2*p2_k_1_SMCKF*A2'+G2*Q2*G2';
        % 对应于各模型进行DMCKF
        t=1;
        sigma=1.1;
      %  sigma=1;
        epsilo=0.01;
        zk=zA(:,k);
        
        Bpk_1 = chol(p1_k_SMCKF)';
        Brk_1 = chol(R)';
        Bkk_1 = blkdiag(Bpk_1,Brk_1);
        dik_1 = inv(Bkk_1)*[x1_pk1_SMCKF;zk];
        wik_1 = [inv(Bpk_1);inv(Brk_1)* H];
        
        xkk_1 = x1_pk1_SMCKF;
        while(t<1000)
            temp=xkk_1;
            t=t+1;
            ee_1 = dik_1 - wik_1 * xkk_1;
            %% 核函数
           % Gee_1 = exp(-ee_1.^2/2/sigma^2);
           % Gee_1 =(dik_1.*(wik_1 * xkk_1)+10).^2; 
             %Gee_1 =dik_1.*(wik_1 * xkk_1);
           %  Gee_1 =tanh(0.2*dik_1.*(wik_1 * xkk_1)+5);
             Gee_1 = 0.9*exp(-ee_1.^2/2/sigma^2)+0.1*tanh(0.2*dik_1.*(wik_1 * xkk_1)+5);
             %Gee_1=0.7*tanh(0.2*dik_1.*(wik_1 * xkk_1)+5)+0.3*dik_1.*(wik_1 * xkk_1);
             %Gee_1 = 0.5*exp(-ee_1.^2/2/sigma^2)+0.5*exp(-ee_1.^2/2/3^2);
             
            Cx_1 = diag(Gee_1(1:4));
            Cy_1 = diag(Gee_1(4+1:6));
%             Pke_hat_1 = inv(Bpk_1)' * (Cx_1) * inv(Bpk_1);
%             R_hat_1 = inv(Brk_1)' * (Cy_1) * inv(Brk_1);
%             Gk_MCC_1 = inv(H'*R_hat_1*H+Pke_hat_1)*H'*R_hat_1;
            
        Pke_hat_1 = Bpk_1 * inv(Cx_1) * Bpk_1';
        R_hat_1 = Brk_1 * inv(Cy_1) * Brk_1';
        Gk_MCC_1 = Pke_hat_1 * H' * inv(H * Pke_hat_1 * H'+R_hat_1);
            
            xkk_1 = x1_pk1_SMCKF + Gk_MCC_1*(zk-H *x1_pk1_SMCKF );
            
            xe_SMCKF_1(:,k) = xkk_1;
            
            if(abs(xe_SMCKF_1(:,k)-temp)/abs(temp)<=epsilo)
                break;
            else
                t=t-1;
            end
        end
        %% 2
        Bpk_2 = chol(p2_k_SMCKF)';
        Brk_2 = chol(R)';
        Bkk_2 = blkdiag(Bpk_2,Brk_2);
        dik_2 = inv(Bkk_2)*[x2_pk1_SMCKF;zk];
        wik_2 = [inv(Bpk_2);inv(Brk_2)* H];
        
        xkk_2 = x2_pk1_SMCKF;
        while(t<1000)
            temp=xkk_2;
            t=t+1;
            ee_2 = dik_2 - wik_2 * xkk_2;
           % Gee_2 = exp(-ee_2.^2/2/sigma^2)
           %Gee_2 =(dik_2.*(wik_2 * xkk_2)+10).^2
            %Gee_2 =dik_2.*(wik_2 * xkk_2)
            %Gee_2 =tanh(0.2*dik_2.*(wik_2 * xkk_2)+5)
           Gee_2 = 0.9*exp(-ee_2.^2/2/sigma^2)+0.1*tanh(0.2*dik_2.*(wik_2 * xkk_2)+5);
            %Gee_2=0.7*tanh(0.2*dik_2.*(wik_2 * xkk_2)+5)+0.3*dik_2.*(wik_2 * xkk_2)
            %Gee_2 = 0.5*exp(-ee_2.^2/2/sigma^2)+0.5*exp(-ee_2.^2/2/3^2);
            Cx_2 = diag(Gee_2(1:4));
            Cy_2 = diag(Gee_2(4+1:6));
%             Pke_hat_2 = inv(Bpk_2)' * (Cx_2) * inv(Bpk_2);
%             R_hat_2 = inv(Brk_2)' * (Cy_2) *inv(Brk_2);
%             Gk_MCC_2 =inv(H'*R_hat_2*H+Pke_hat_2)*H'*R_hat_2;
            
        Pke_hat_2 = Bpk_2 * inv(Cx_2) * Bpk_2';
        R_hat_2 = Brk_2 * inv(Cy_2) * Brk_2';
        Gk_MCC_2 = Pke_hat_2 * H' * inv(H * Pke_hat_2 * H'+R_hat_2);
            
        
            xkk_2 = x2_pk1_SMCKF + Gk_MCC_2*(zk-H *x2_pk1_SMCKF );
            
            xe_SMCKF_2(:,k) = xkk_2;
            
            if(abs(xe_SMCKF_2(:,k)-temp)/abs(temp)<=epsilo)
                break;
            else
                t=t-1;
            end
        end
     
        %预测残差及协方差计算
        v1_SMCKF=zk-H*x1_pk1_SMCKF;
        v2_SMCKF=zk-H*x2_pk1_SMCKF;
        Sv1_SMCKF=H*p1_k_SMCKF*H'+R;
        Sv2_SMCKF=H*p2_k_SMCKF*H'+R;
        like1_SMCKF=det(2*pi*Sv1_SMCKF)^(-0.5)*exp(-0.5*v1_SMCKF'*inv(Sv1_SMCKF)*v1_SMCKF);
        like2_SMCKF=det(2*pi*Sv2_SMCKF)^(-0.5)*exp(-0.5*v2_SMCKF'*inv(Sv2_SMCKF)*v2_SMCKF);
        %滤波更新
        K1_SMCKF=Gk_MCC_1;K2_SMCKF=Gk_MCC_2;
        xk1_SMCKF=xe_SMCKF_1(:,k);xk2_SMCKF=xe_SMCKF_2(:,k);
        P1_SMCKF=(eye(4) - Gk_MCC_1 * H) * p1_k_SMCKF * (eye(4) - Gk_MCC_1 * H)'+ Gk_MCC_1 * R * Gk_MCC_1';
        P2_SMCKF=(eye(4) - Gk_MCC_2 * H) * p2_k_SMCKF * (eye(4) - Gk_MCC_2 * H)'+ Gk_MCC_2 * R * Gk_MCC_2';
        %模型概率更新
        C_SMCKF=like1_SMCKF*c1_SMCKF+like2_SMCKF*c2_SMCKF;
        u1_SMCKF=like1_SMCKF*c1_SMCKF/C_SMCKF;
        u2_SMCKF=like2_SMCKF*c2_SMCKF/C_SMCKF;
        %估计融合
        xk_SMCKF=xk1_SMCKF*u1_SMCKF+xk2_SMCKF*u2_SMCKF;
        %迭代
        x1_k_1_SMCKF=xk1_SMCKF;x2_k_1_SMCKF=xk2_SMCKF;
        X_imm_SMCKF(:,k)=xk_SMCKF;
        um1_SMCKF(k)=u1_SMCKF;um2_SMCKF(k)=u2_SMCKF;
    end
    
    ex_mckf_imm(mn,:)=X_imm_MCKF(1,:)-xA(1,:);
    ey_mckf_imm(mn,:)=X_imm_MCKF(3,:)-xA(3,:);
    ex_smckf_imm(mn,:)=X_imm_SMCKF(1,:)-xA(1,:);
    ey_smckf_imm(mn,:)=X_imm_SMCKF(3,:)-xA(3,:);
    
    m_x_xA(mn,:)=xA(1,:);
    m_y_xA(mn,:)=xA(3,:);
    m_x_zA(mn,:)=zA(1,:);
    m_y_zA(mn,:)=zA(2,:);
    
    
    UM1_MCKF(mn,:)=um1_MCKF;
    UM2_MCKF(mn,:)=um2_MCKF;
    UM1_SMCKF(mn,:)=um1_SMCKF;
    UM2_SMCKF(mn,:)=um2_SMCKF;
end

mex_kf_cv=mean(ex_kf_cv);mey_kf_cv=mean(ey_kf_cv);%CV
mex_kf_ct=mean(ex_kf_ct);mey_kf_ct=mean(ey_kf_ct);%CT
mex_mckf_cv=mean(ex_mckf_cv);mey_mckf_cv=mean(ey_mckf_cv);%MCKF_CV
mex_mckf_ct=mean(ex_mckf_ct);mey_mckf_ct=mean(ey_mckf_ct);%MCKF_CT
mex_imm=mean(ex_imm);mey_imm=mean(ey_imm);%IMM
mex_mckf_imm=mean(ex_mckf_imm);mey_mckf_imm=mean(ey_mckf_imm);%IMM-mckf
mex_smckf_imm=mean(ex_smckf_imm);mey_smckf_imm=mean(ey_smckf_imm);%imm-smckf
% mm_x_xA=mean(m_x_xA);mm_y_xA=mean(m_y_xA);
% mm_x_zA=mean(m_x_zA);mm_y_zA=mean(m_y_zA);

%% 绘图
t=1:Len;

%% 绘制各算法跟踪效果
figure(1)
plot(xA(1,t),xA(3,t),'k' ,zA(1,t),zA(2,t),'m',xA(1,t)+mex_imm(t),xA(3,t)+mey_imm(t),...
    'g--',xA(1,t)+mex_mckf_cv(t),xA(3,t)+mey_mckf_cv(t),'c--',xA(1,t)+mex_mckf_imm(t),xA(3,t)+mey_mckf_imm(t),'b--',xA(1,t)+mex_smckf_imm(t),xA(3,t)+mey_smckf_imm(t),'r--','linewidth',1);

% plot(xA(1,t),xA(3,t),'k' ,zA(1,t),zA(2,t),'m',X_imm(1,t),X_imm(3,t),...
%     'r--',X_imm_MCKF(1,t),X_imm_MCKF(3,t),'g--',X_MCKF_CV(1,t),X_MCKF_CV(3,t),'b--',X_imm_SMCKF(1,t),X_imm_SMCKF(3,t),'c--','linewidth',1);

legend('true', 'measurment','IMM-KF','MCKF-CV','IMM-MCKF','IMM-IMCKF');
title('Tracking results')
xlabel('x/(m)'),ylabel('y/(m)');

%% 误差图
figure(2)
subplot(2,1,1)
plot(t,mex_kf_cv(t),'b.-',t,mex_kf_ct(t),'g.-',t,mex_imm(t),'r.-',t,mex_mckf_imm(t),'k.-',t,mex_smckf_imm(t),'c.-');
%plot(t,20*log10(mex1(t)/max(mex1)),'b',t,20*log10(mex2(t)/max(mex2)),'m',t,20*log10(mex3(t)/max(mex3)),'r',t,20*log10(mex3_MCKF(t)/max(mex3_MCKF)),'y');
title('X方向滤波误差')
legend('CV模型滤波','CT模型滤波','IMM滤波','IMM-MCKF滤波','IMM-SMCKF滤波');
subplot(2,1,2)
plot(t,mey_kf_cv(t),'b.-',t,mey_kf_ct(t),'g.-',t,mey_imm(t),'r.-',t,mey_mckf_imm(t),'k.-',t,mey_smckf_imm(t),'c.-');
%plot(t,20*log10(mey1(t)/max(mey1)),'b',t,20*log10(mey2(t)/max(mey2)),'m',t,20*log10(mey3(t)/max(mey3)),'r',t,20*log10(mey3_MCKF(t)/max(mey3_MCKF)),'y');
legend('CV模型滤波','CT模型滤波','IMM滤波','IMM-MCKF滤波','IMM-SMCKF滤波');
title('Y方向滤波误差')
xlabel('t(s)'),ylabel('位置误差(m)');

%%  计算RMSE-均方根误差
EX_kf_cv=sqrt(sum(ex_kf_cv.^2)/MC);
EX_kf_ct=sqrt(sum(ex_kf_ct.^2)/MC);
EX_imm=sqrt(sum(ex_imm.^2)/MC);
EX_mckf_imm=sqrt(sum(ex_mckf_imm.^2)/MC);
EX_smckf_imm=sqrt(sum(ex_smckf_imm.^2)/MC);
EX_mckf_cv=sqrt(sum(ex_mckf_cv.^2)/MC);
EX_mckf_ct=sqrt(sum(ex_mckf_ct.^2)/MC);

EY_kf_cv=sqrt(sum(ey_kf_cv.^2)/MC);
EY_kf_ct=sqrt(sum(ey_kf_ct.^2)/MC);
EY_imm=sqrt(sum(ey_imm.^2)/MC);
EY_mckf_imm=sqrt(sum(ey_mckf_imm.^2)/MC);
EY_smckf_imm=sqrt(sum(ey_smckf_imm.^2)/MC);
EY_mckf_cv=sqrt(sum(ey_mckf_cv.^2)/MC);
EY_mckf_ct=sqrt(sum(ey_mckf_ct.^2)/MC);


figure(3)
subplot(1,2,1)
%plot( t,EX1(t),'b.-',t,EX2(t),'k.-',t,EX3(t),'r.-',t,EX3_MCKF(t),'g.-');
%plot( t,EX1(t),'b.-',t,EX2(t),'k.-',t,EX3(t),'r.-',t,EX3_MCKF(t),'g.-',t,EX5(t),'y.-',t,EX6(t),'m.-');
%plot(t,EX3(t),'r.-',t,EX3_MCKF(t),'g.-',t,EX5(t),'b.-');
plot(t,20*log10(EX_imm(t)),'g.-',t,20*log10(EX_mckf_cv(t)),'c.-',t,20*log10(EX_mckf_imm(t)),'b.-',t,20*log10(EX_smckf_imm(t)),'r.-');
title('X-direction location RMSE')
xlabel('time/(s)'),ylabel('RMSE/(db)');
%xlabel('time/(s)'),ylabel('RMSE/(m)');
%legend('CV-KF','CT-KF','IMM','DMC-IMM','CV-MCKF','CT-MCKF');
legend('IMM-KF','MCKF-CV','IMM-MCKF','IMM-IMCKF');
mean_x_cv=mean(EX_kf_cv(1,400:end));mean_x_ct=mean(EX_kf_ct(1,400:end));mean_x_imm=mean(EX_imm(1,400:end));mean_x_IMM_MCKF=mean(EX_mckf_imm(1,400:end));mean_x_IMM_SMCKF=mean(EX_smckf_imm(1,400:end));
mean_x_MCKF=mean(EX_mckf_cv(1,400:end));
%figure(4)
subplot(1,2,2)
%plot(t,EY1(t),'b.-',t,EY2(t),'k.-',t,EY3(t),'r.-',t,EY3_MCKF(t),'g.-',t,EY5(t),'y.-',t,EY6(t),'m.-' );
%legend('CV-KF','CT-KF','IMM','DMC-IMM','CV-MCKF','CT-MCKF');
plot(t,20*log10(EY_imm(t)),'g.-',t,20*log10(EY_mckf_cv(t)),'c.-',t,20*log10(EY_mckf_imm(t)),'b.-',t,20*log10(EY_smckf_imm(t)),'r.-');
legend('IMM-KF','MCKF-CV','IMM-MCKF','IMM-IMCKF');
title('Y-direction location RMSE')
%xlabel('time/(s)'),ylabel('RMSE/(m)');
xlabel('time/(s)'),ylabel('RMSE/(db)');
mean_y_cv=mean(EY_kf_cv(1,400:end));mean_y_ct=mean(EY_kf_ct(1,400:end));mean_y_imm=mean(EY_imm(1,400:end));mean_y_IMM_MCKF=mean(EY_mckf_imm(1,400:end));mean_y_IMM_SMCKF=mean(EY_smckf_imm(1,400:end));
mean_y_MCKF=mean(EY_mckf_cv(1,400:end));


Um1=mean(UM1);Um2=mean(UM2);
Um1_MCKF=mean(UM1_MCKF);Um2_MCKF=mean(UM2_MCKF);%各模型概率
Um1_SMCKF=mean(UM1_SMCKF);Um2_SMCKF=mean(UM2_SMCKF);
figure(4)
plot(t,Um1(t),'b-',t,Um2(t),'m-',t,Um1_MCKF(t),'k-',t,Um2_MCKF(t),'y-',t,Um1_SMCKF(t),'c-',t,Um2_SMCKF(t),'g-' );
legend('IMM-CV模型概率','IMM-CT模型概率','IMM-MCKF-CV模型概率','IMM-MCKF-CT模型概率','IMM-SMCKF-CV模型概率','IMM-SMCKF-CT模型概率');
title('CV、CT模型概率变化')
toc
