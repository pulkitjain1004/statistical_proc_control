clear all
filename = 'Data.xls';
sheet = 1;
xlrange = 'A1:HA552';
Datam = xlsread(filename, sheet, xlrange, 'basic'); 
[m n] = size(Datam);
Datam = Datam';

%calculate mean
xbar = sum(Datam,2)/m;
%disp(xbar)

%calcualte covariance matrix S
S2 = zeros(n);
for i=1:m
    S2 = S2 + (Datam(:,i) - xbar)*(Datam(:,i) - xbar)';
end
S = S2/(m-1);
diagS = diag(S);

[E Lambda] = eig(S);
E = E(:, n:-1:1);
Lambda = Lambda(n:-1:1, n:-1:1);

%calulate correlation matrix rho
rho = zeros(209);
for i = 1:209
   for j=1:209
      rho(i,j) = S(i,j)/sqrt(S(i,i)*S(j,j));
   end
end
  
%calculate eigen values of correltaion matrix rho
[Ep Lambda_rho ] = eig(rho);
 Ep = Ep(:, n:-1:1);
 Lambda_rho = Lambda_rho(n:-1:1, n:-1:1);
 lambda_rho = diag(Lambda_rho);
 
%create scree plot
%figure(1)
plot(lambda_rho(1:10));
grid on
grid minor
title('Scree plot for Correlation Matrix')
xlabel('Index Value')
ylabel('Lambda_rho')

%create pareto plot
%figure(2)
pareto(lambda_rho);
grid on
grid minor
title('Pareto plot for Correlation Matrix')
xlabel('Index Value')
ylabel('Lambda_rho')


format long
sample_size = 552;
 MDL_rho = zeros(208,1);
 r=1;
 for h=0:208
     a = mean(lambda_rho(h+1:n));
     g = geomean(lambda_rho(h+1:n));
     MDL_rho(r) = sample_size*(n-h)*log(a/g)+h*(2*n-h)*log(sample_size)/2;
     r=r+1;
 end

% figure(3)
plot(MDL_rho);
grid on
grid minor
title('MDL plot for Correlation Matrix')
xlabel('Index Value')
ylabel('MDL')

 %calculate V^1/2 matrix
 Vhalf = zeros(209);
 for w=1:209
     Vhalf(w,w)= 1/sqrt(S(w,w));
 end
 
 %PCA for correlation matrix where the ith PC is yi=ei transpose
PC1p = zeros(1,m);
PC2p = zeros(1,m);
PC3p = zeros(1,m);
PC4p = zeros(1,m);
PC5p = zeros(1,m);
 for w = 1:m
     PC1p(w) = Ep(:,1)'*Vhalf*(Datam(:,w)-xbar);
     PC2p(w) = Ep(:,2)'*Vhalf*(Datam(:,w)-xbar);
     PC3p(w) = Ep(:,3)'*Vhalf*(Datam(:,w)-xbar);
     PC4p(w) = Ep(:,4)'*Vhalf*(Datam(:,w)-xbar);
     PC5p(w) = Ep(:,5)'*Vhalf*(Datam(:,w)-xbar);
 end

UCL1p = 3*sqrt(Lambda_rho(1,1));
UCL2p = 3*sqrt(Lambda_rho(2,2));
UCL3p = 3*sqrt(Lambda_rho(3,3));
UCL4p = 3*sqrt(Lambda_rho(4,4));
UCL5p = 3*sqrt(Lambda_rho(5,5));

LCL1p = -3*sqrt(Lambda_rho(1,1));
LCL2p = -3*sqrt(Lambda_rho(2,2));
LCL3p = -3*sqrt(Lambda_rho(3,3));
LCL4p = -3*sqrt(Lambda_rho(4,4));
LCL5p = -3*sqrt(Lambda_rho(5,5));

%plot individual x-control charts
figure(4)
plot(1:m,PC1p)
hold on
title('PC1p')
xlabel('Index')
xL = get(gca,'XLim');
line(xL,[UCL1p UCL1p],'Color','r');
line(xL,[LCL1p LCL1p],'Color','r');
hold off

figure(5)
plot(1:m,PC2p)
hold on
title('PC2p')
xL = get(gca,'XLim');
line(xL,[UCL2p UCL2p],'Color','r');
line(xL,[LCL2p LCL2p],'Color','r');
hold off

figure(6)
plot(1:m,PC3p)
hold on
title('PC3p')
xL = get(gca,'XLim');
line(xL,[UCL3p UCL3p],'Color','r');
line(xL,[LCL3p LCL3p],'Color','r');
hold off

figure(7)
plot(1:m,PC4p)
hold on
title('PC4p')
xL = get(gca,'XLim');
line(xL,[UCL4p UCL4p],'Color','r');
line(xL,[LCL4p LCL4p],'Color','r');
hold off

figure(8)
plot(1:m,PC5p)
hold on
title('PC5p')
xL = get(gca,'XLim');
line(xL,[UCL5p UCL5p],'Color','r');
line(xL,[LCL5p LCL5p],'Color','r');
hold off

A = [PC1p;PC2p;PC3p;PC4p;PC5p];
%plot along with Tsquare & mCUSUM

[JJ J] = size(A);  
Abar = sum(A,2)/J;
S3 = zeros(5);
for i=1:J
    S3 = S3 + (A(:,i) - Abar)*(A(:,i) - Abar)';
end
A_cov = S3/(J-1);
A_cov_inv = inv(A_cov);

UCL = chi2inv(0.995,5);
    for j = 1:J
        T(j) = (A(:,j) - Abar)'*A_cov_inv*(A(:,j) - Abar);
    end
       
    figure(9)
    plot(T(1:J))
    hold on
    title('T^2 Chart for Correlation for alpha = 0.005, entire range of observation')
    xlabel('Data Index')
    ylabel('T2 Statisitc')
    xL = get(gca,'XLim');
    line(xL,[UCL UCL],'Color','r');
    text(m+0.5, UCL, 'UCL');
    hold off

% mCUSUM entire range   

UCL_now = 6.6;  %p=5 ARL = 200 interpolate
C_i = zeros(5, J);
    MC_i = zeros(1, J);
    n_i = zeros(1, J);
    x_pro = A(:,1:J) - repmat(Abar, 1, J);
    
    k = 0.5 * 3;
        for i = 1:J
        if i == 1
            n_i(i) = 1;
        else
            if MC_i(i-1) > 0
                n_i(i) = n_i(i-1) + 1;
            else
                n_i(i) = 1;
            end
        end
        j = i - n_i(i) + 1;
        C_i(:, i) = sum(x_pro(:, j:i), 2);
        MC_i(i) = max(0, sqrt(C_i(:, i)' * A_cov_inv * C_i(:, i)) - k * n_i(i));
        end
  
       
    figure(10)
    plot(MC_i(1:J))
    hold on
    title('mCUSUM Chart for 3 sigma mean shift, Correlation for entire range')
    xlabel('Data Index')
    ylabel('T2 Statisitc')
    xL = get(gca,'XLim');
    line(xL,[UCL_now UCL_now],'Color','r');
    text(m+0.5, UCL_now, 'UCL');
    hold off
    
    
B = horzcat(A(:,52:269),A(:,271:445),A(:,465:529)); %feed in manually
[JJ J] = size(B);

warm = Datam(:,1:51);
Warmup = sum(warm,2)/51;
seg2 = Datam(:,270);
Seg2 = seg2;
seg3 = Datam(:,530:552);
Seg3 = sum(seg3,2)/(552-530+1);
seg4 = Datam(:,446:464);
Seg4 = sum(seg4,2)/(464-446+1);
b_temp = horzcat(Datam(:,52:269),Datam(:,271:445),Datam(:,465:529));  %,Datam(:,538:552));
B_temp = sum(b_temp,2)/J;

figure(1004)
plot(1:209,Warmup,'color','r');
hold on
title('Original Signal for 4 segments')
xlabel('Input Variables')
ylabel('Output Statisitc')
plot(1:209,Seg2,'color','c')
plot(1:209,Seg3,'color','g')
plot(1:209,Seg4,'color','g')
plot(1:209,B_temp,'color','b')
hold off

%calculate Hotelling T^2 for covariance matrix PCA reduction data

[JJ J] = size(B);
J
bbar_cov = sum(B,2)/J;
S3 = zeros(5);
for i=1:J
    S3 = S3 + (B(:,i) - bbar_cov)*(B(:,i) - bbar_cov)';
end
B_cov = S3/(J-1);
B_cov_inv = inv(B_cov);

UCL = chi2inv(0.995,5);

check = 0;
ii=1;
while(check~=bbar_cov)
    Tnew = zeros(5,J);
    check = bbar_cov;
    q=1;
    for j = 1:J
        T(j) = (B(:,j) - bbar_cov)'*B_cov_inv*(B(:,j) - bbar_cov);
        if(T(j) < UCL)
         Tnew(:,q)=B(:,j);
         q=q+1;
        end
    end
    
    ii=ii+1;
    
    figure(230+ii)
    plot(T(1:J))
    hold on
    title('T^2 Chart for alpha = 0.005')
    xlabel('Data Index')
    ylabel('T2 Statisitc')
    xL = get(gca,'XLim');
    line(xL,[UCL UCL],'Color','r');
    text(m+0.5, UCL, 'UCL');
    hold off 

    q=q-1;
    xbarnew = sum(Tnew,2)/q;
    bbar_cov = xbarnew;
    
    J = q
    S2new = zeros(5);       %recalcualte covariance matrix S
        for i=1:J
            S2new = S2new + (Tnew(:,i) - bbar_cov)*(Tnew(:,i) - bbar_cov)';
        end
    Snew = S2new/(J-1);
    B_cov_inv = inv(Snew);
    B = Tnew;

end

%m-CUSUM

D = B(:,1:J);      % new matric B   
m =J;
Dbar = sum(D,2)/m;
check2 = 0;
qqq=1;
while(check2~=Dbar)
G = zeros(5,552);
qq=1;
check2 = Dbar;
DS = zeros(5);       %recalcualte covariance matrix S
    for i=1:m
        DS = DS + (D(:,i) - Dbar)*(D(:,i) - Dbar)';
    end
DS = DS/(m-1);

% disp(Snew)
DS_inv = inv(DS);
UCL_now = 6.6;  %p=5 ARL = 200
C_i = zeros(5, m);
    MC_i = zeros(1, m);
    n_i = zeros(1, m);
    x_pro = D(:,1:m) - repmat(Dbar, 1, m);
    
    k = 0.5 * 3;
        for i = 1:m
        if i == 1
            n_i(i) = 1;
        else
            if MC_i(i-1) > 0
                n_i(i) = n_i(i-1) + 1;
            else
                n_i(i) = 1;
            end
        end
        j = i - n_i(i) + 1;
        C_i(:, i) = sum(x_pro(:, j:i), 2);
        MC_i(i) = max(0, sqrt(C_i(:, i)' * DS_inv * C_i(:, i)) - k * n_i(i));
        end
  
    for i = 1:m
        if MC_i(i)<UCL_now
            G(:,qq)=D(:,i);
            qq = qq+1;
        end
    end

    qqq=qqq+1;
    
    figure(700+qqq)
    plot(MC_i(1:m))
    hold on
    title('mCUSUM Chart for 3 sigma mean shift')
    xlabel('Data Index')
    ylabel('T2 Statisitc')
    xL = get(gca,'XLim');
    line(xL,[UCL_now UCL_now],'Color','r');
    text(m+0.5, UCL_now, 'UCL');
    hold off
    
    qq=qq-1
    D = G;
    m=qq;
    Dbar = sum(D,2)/m;
end 
    
    
%Round 2 calculate Hotelling T^2 for covariance matrix PCA reduction

F = D(:,1:m);

[LL L] = size(F);  
Fbar = sum(F,2)/L;
F3 = zeros(5);
for i=1:L
    F3 = F3 + (F(:,i) - Fbar)*(F(:,i) - Fbar)';
end
F_cov = F3/(L-1);
F_cov_inv = inv(F_cov);


UCL = chi2inv(0.995,5);

check = 0;
ii=20;
while(check~=Fbar)
    Fnew = zeros(5,J);
    check = Fbar;
    q=1;
    for j = 1:L
        T_F(j) = (F(:,j) - Fbar)'*F_cov_inv*(F(:,j) - Fbar);
        if(T_F(j) < UCL)
         Fnew(:,q)=F(:,j);
         q=q+1;
        end
    end
    
    ii=ii+1;
    
    figure(230+ii)
    plot(T_F(1:L))
    hold on
    title('Round 2 - T^2 Chart for alpha = 0.005')
    xlabel('Data Index')
    ylabel('T2 Statisitc')
    xL = get(gca,'XLim');
    line(xL,[UCL UCL],'Color','r');
    text(m+0.5, UCL, 'UCL');
    hold off 

    q=q-1;
    Fbarnew = sum(Fnew,2)/q;
    Fbar = Fbarnew;
    
    L = q
    F2new = zeros(5);       %recalcualte covariance matrix S
        for i=1:L
            F2new = F2new + (Fnew(:,i) - Fbar)*(Fnew(:,i) - Fbar)';
        end
    F2new = F2new/(L-1);
    % disp(Snew)
    F_cov_inv = inv(F2new);
    F = Fnew;

end

%m-CUSUM round 2

H = F(:,1:L);      % new matric B   
m =L;
Hbar = sum(H,2)/m;
check4 =0;
qqq=10;

while(check4~=Hbar)
G2 = zeros(5,552);
qq2=1;
check4 = Hbar;
HS = zeros(5);       %recalcualte covariance matrix S
    for i=1:m
        HS = HS + (H(:,i) - Hbar)*(H(:,i) - Hbar)';
    end
HS = HS/(m-1);

HS_inv = inv(HS);
UCL_now = 6.6;  %p=5 ARL = 200
C_i = zeros(5, m);
    MC_i = zeros(1, m);
    n_i = zeros(1, m);
    x_pro = H(:,1:m) - repmat(Hbar, 1, m);
    
    k = 0.5 * 3;
        for i = 1:m
        if i == 1
            n_i(i) = 1;
        else
            if MC_i(i-1) > 0
                n_i(i) = n_i(i-1) + 1;
            else
                n_i(i) = 1;
            end
        end
        j = i - n_i(i) + 1;
        C_i(:, i) = sum(x_pro(:, j:i), 2);
        MC_i(i) = max(0, sqrt(C_i(:, i)' * HS_inv * C_i(:, i)) - k * n_i(i));
        end
  
    for i = 1:m
        if MC_i(i)<UCL_now
            G2(:,qq2)=H(:,i);
            qq2 = qq2+1;
        end
    end

    qqq = qqq+1;
    figure(800+qqq)
    plot(MC_i(1:m))
    hold on
    title('Round 2 mCUSUM Chart for 3 sigma mean shift')
    xlabel('Data Index')
    ylabel('MC Statisitc')
    xL = get(gca,'XLim');
    line(xL,[UCL_now UCL_now],'Color','r');
    text(m+0.5, UCL_now, 'UCL');
    hold off
    
    qq2=qq2-1
    H=G2;
    m=qq2;
    Hbar = sum(H,2)/m;
end

Final = H(:,1:m);
[LL L] = size(Final);  
Finalbar = sum(Final,2)/L
F3 = zeros(5);
for i=1:L
    F3 = F3 + (Final(:,i) - Finalbar)*(Final(:,i) - Fbar)';
end
Final_cov = F3/(L-1)
F_cov_inv = inv(F_cov); %For PC1, PC3 & PC3