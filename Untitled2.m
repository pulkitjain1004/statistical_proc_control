clear all
filename = 'Data.xls';
sheet = 1;
xlrange = 'A1:HA552';
Datam = xlsread(filename, sheet, xlrange, 'basic'); 
[m n] = size(Datam);
Datam = Datam';

%calculate mean
xbar = sum(Datam,2)/m;

%calcualte covariance matrix S
S2 = zeros(n);
for i=1:m
    S2 = S2 + (Datam(:,i) - xbar)*(Datam(:,i) - xbar)';
end
S = S2/(m-1);
diagS = diag(S);

%plot the original variance to see how much it differs between values
figure(1)
plot(diagS , 'Linewidth', 1)
title('Plot for Individual Variance Of Original Components','fontweight','bold','fontsize',14)
xlabel('Component Index','fontweight','bold')
ylabel('Variance','fontweight','bold')

%calculate eigenvectors and eigenvalues of S
[E Lambda] = eig(S);
%reaarange in descending order
E = E(:, n:-1:1);
Lambda = Lambda(n:-1:1, n:-1:1);

%calulate MDL for lambda when n=1
lambda = diag(Lambda);
format long
sample_size = 552;
MDL = zeros(208,1);
r=1;
for h=0:208
    a = mean(lambda(h+1:n));
    g = geomean(lambda(h+1:n));
    MDL(r) = sample_size*(n-h)*log(a/g)+h*(2*n-h)*log(sample_size)/2;
    r=r+1;
end

%plot MDL
 figure(2)
 plot(0:208, MDL, 'Linewidth', 1);
 %grid on
 %grid minor
 title('MDL Plot for Covariance Matrix','fontweight','bold','fontsize',14)
 xlabel('Principal Component Index','fontweight','bold')
 ylabel('MDL','fontweight','bold')

%create scree plot
figure(3)
plot(lambda(1:10), 'Linewidth', 1);
grid on
%grid minor
title('Scree plot for Covariance Matrix','fontweight','bold','fontsize',14)
xlabel('Principal Component Index','fontweight','bold')
ylabel('Variance / Lambda','fontweight','bold')


%create pareto plot
figure(4)
pareto(lambda)
grid on
%grid minor
title('Pareto plot for Covariance Matrix','fontweight','bold','fontsize',14)
xlabel('Principal Component Index','fontweight','bold')
ylabel('Variance / Lambda','fontweight','bold')


%PCA for covariance matrix where the ith PC is yi=ei transpose*(x-xbar)
PC1 = zeros(1,m);
PC2 = zeros(1,m);
PC3 = zeros(1,m);
for i=1:m
    PC1(i) = E(:,1)'*(Datam(:,i)-xbar);
    PC2(i) = E(:,2)'*(Datam(:,i)-xbar);
    PC3(i) = E(:,3)'*(Datam(:,i)-xbar);
end
UCL1 = 3*sqrt(Lambda(1,1));
UCL2 = 3*sqrt(Lambda(2,2));
UCL3 = 3*sqrt(Lambda(3,3));

LCL1 = -3*sqrt(Lambda(1,1));
LCL2 = -3*sqrt(Lambda(2,2));
LCL3 = -3*sqrt(Lambda(3,3));


%plot individual PC charts

figure(5)
plot(1:m,PC1, 'Linewidth', 1);
hold on
title('First Principal Component','fontweight','bold','fontsize',14)
xlabel('Index')
ylabel('PC 1','fontweight','bold')
xL = get(gca,'XLim');
line(xL,[UCL1 UCL1],'Color','r');
line(xL,[LCL1 LCL1],'Color','r');
xlabel('Data Index','fontweight','bold')
text(553, UCL1, 'UCL');
text(553, LCL1, 'LCL');
hold off

figure(6)
plot(1:m,PC2, 'Linewidth', 1);
hold on
title('Second Principal Component','fontweight','bold','fontsize',14)
xL = get(gca,'XLim');
line(xL,[UCL2 UCL2],'Color','r');
line(xL,[LCL2 LCL2],'Color','r');
xlabel('Data Index','fontweight','bold')
ylabel('PC 2','fontweight','bold')
text(553, UCL2, 'UCL');
text(553, LCL2, 'LCL');
hold off

figure(7)
plot(1:m,PC3, 'Linewidth', 1);
hold on
title('Third Principal Component','fontweight','bold','fontsize',14)
xL = get(gca,'XLim');
line(xL,[UCL3 UCL3],'Color','r');
line(xL,[LCL3 LCL3],'Color','r');
xlabel('Data Index','fontweight','bold')
ylabel('PC 3','fontweight','bold')
text(553, UCL3, 'UCL');
text(553, LCL3, 'LCL');
hold off

%plot along with Tsquare & mCUSUM

A = [PC1;PC2;PC3];
B = horzcat(A(:,49:443),A(:,458:529)); %,A(:,538:552)); %feed in manually

warm = Datam(:,1:48);
Warmup = sum(warm,2)/48;
seg2 = Datam(:,444:457);
Seg2 = sum(seg2,2)/(457-444+1);
seg3 = Datam(:,530:552);
Seg3 = sum(seg3,2)/(552-530+1);
b_temp = horzcat(Datam(:,49:443),Datam(:,458:529));  %,Datam(:,538:552));
B_temp = sum(b_temp,2)/467;

figure(10)
plot(1:209,Warmup,'-xr');
hold on
title('Original Signal for 4 segments','fontweight','bold','fontsize',14)
xlabel('Input Variables','fontweight','bold')
ylabel('Output Statisitc','fontweight','bold')
plot(1:209,Seg2,'-xc')
plot(1:209,Seg3,'-xg')
plot(1:209,B_temp,'-ob')
hold off

%calculate Hotelling T^2 for covariance matrix PCA reduction data

[JJ J] = size(B);
J
bbar_cov = sum(B,2)/J;
S3 = zeros(3);
for i=1:J
    S3 = S3 + (B(:,i) - bbar_cov)*(B(:,i) - bbar_cov)';
end
B_cov = S3/(J-1);
B_cov_inv = inv(B_cov);

% Round 1, calculate T^2

UCL = chi2inv(0.995,3);

check = 0;
ii=1;
while(check~=bbar_cov)
    Tnew = zeros(3,J);
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
    
    figure(98+ii)
    plot(T(1:J), 'Linewidth', 1)
    hold on
    title('Round-1: T^2 Chart for alpha = 0.005, Iteration no.','fontweight','bold','fontsize',14)
    xlabel('Data Index','fontweight','bold')
    ylabel('T2 Statisitc','fontweight','bold')
    xL = get(gca,'XLim');
    line(xL,[UCL UCL],'Color','r');
    text(225, UCL, 'UCL');
    hold off 

    q=q-1;
    xbarnew = sum(Tnew,2)/q;
    bbar_cov = xbarnew;
    
    J = q
    S2new = zeros(3);       %recalcualte covariance matrix S
        for i=1:J
            S2new = S2new + (Tnew(:,i) - bbar_cov)*(Tnew(:,i) - bbar_cov)';
        end
    Snew = S2new/(J-1);
    
    B_cov_inv = inv(Snew);
    B = Tnew;

end

    
%round 1, m-CUSUM

D = B(:,1:J);      % new matric B   
m =J;
Dbar = sum(D,2)/m;
check2 = 0;
qqq=1;
while(check2~=Dbar)
%disp(C);
%plot(D_T);
G = zeros(3,552);
qq=1;
check2 = Dbar;
DS = zeros(3);       %recalcualte covariance matrix S
    for i=1:m
        DS = DS + (D(:,i) - Dbar)*(D(:,i) - Dbar)';
    end
DS = DS/(m-1);

% disp(Snew)
DS_inv = inv(DS);
UCL_now = 5.48;  %p=3 ARL = 199.74
C_i = zeros(3, m);
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
    
    figure(198+qqq)
    plot(MC_i(1:m), 'Linewidth', 1)
    hold on
    title('Round-1: mCUSUM Chart for 3 sigma mean shift, Iteration no. ','fontweight','bold','fontsize',14)
    xlabel('Data Index','fontweight','bold')
    ylabel('MC Statisitc','fontweight','bold')
    xL = get(gca,'XLim');
    line(xL,[UCL_now UCL_now],'Color','r');
    text(220, UCL_now, 'UCL');
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
F3 = zeros(3);
for i=1:L
    F3 = F3 + (F(:,i) - Fbar)*(F(:,i) - Fbar)';
end
F_cov = F3/(L-1);
F_cov_inv = inv(F_cov);


UCL = chi2inv(0.995,3);

check = 0;
ii=20;
while(check~=Fbar)
    Fnew = zeros(3,J);
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
    
    figure(279+ii)
    plot(T_F(1:L), 'Linewidth', 1)
    hold on
    title('Round-2: T^2 Chart for alpha = 0.005, Iteration no. ','fontweight','bold','fontsize',14)
    xlabel('Data Index','fontweight','bold')
    ylabel('T2 Statisitc','fontweight','bold')
    xL = get(gca,'XLim');
    line(xL,[UCL UCL],'Color','r');
    text(218, UCL, 'UCL');
    hold off 

    q=q-1;
    Fbarnew = sum(Fnew,2)/q;
    Fbar = Fbarnew;
    
    L = q
    F2new = zeros(3);       %recalcualte covariance matrix S
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
G2 = zeros(3,552);
qq2=1;
check4 = Hbar;
HS = zeros(3);       %recalcualte covariance matrix S
    for i=1:m
        HS = HS + (H(:,i) - Hbar)*(H(:,i) - Hbar)';
    end
HS = HS/(m-1);

HS_inv = inv(HS);
UCL_now = 5.48;  %p=3 ARL = 199.74
C_i = zeros(3, m);
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
    figure(389+qqq)
    plot(MC_i(1:m), 'Linewidth', 1)
    hold on
    title('Round-2: mCUSUM Chart for 3 sigma mean shift, Iteration no. ','fontweight','bold','fontsize',14)
    xlabel('Data Index','fontweight','bold')
    ylabel('MC Statisitc','fontweight','bold')
    xL = get(gca,'XLim');
    line(xL,[UCL_now UCL_now],'Color','r');
    text(218, UCL_now, 'UCL');
    hold off
    
    qq2=qq2-1
    H=G2;
    m=qq2;
    Hbar = sum(H,2)/m;
end

%Round 3 T square

F = H(:,1:m);

[LL L] = size(F);  
Fbar = sum(F,2)/L;
F3 = zeros(3);
for i=1:L
    F3 = F3 + (F(:,i) - Fbar)*(F(:,i) - Fbar)';
end
F_cov = F3/(L-1);
F_cov_inv = inv(F_cov);


UCL = chi2inv(0.995,3);

check = 0;
ii=30;
while(check~=Fbar)
    Fnew = zeros(3,J);
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
    
    figure(469+ii)
    plot(T_F(1:L), 'Linewidth', 1)
    hold on
    title('Round-3: T^2 Chart for alpha = 0.005, Iteration no. ','fontweight','bold','fontsize',14)
    xlabel('Data Index','fontweight','bold')
    ylabel('T2 Statisitc','fontweight','bold')
    xL = get(gca,'XLim');
    line(xL,[UCL UCL],'Color','r');
    text(218, UCL, 'UCL');
    hold off 

    q=q-1;
    Fbarnew = sum(Fnew,2)/q;
    Fbar = Fbarnew;
    
    L = q
    F2new = zeros(3);       %recalcualte covariance matrix S
        for i=1:L
            F2new = F2new + (Fnew(:,i) - Fbar)*(Fnew(:,i) - Fbar)';
        end
    F2new = F2new/(L-1);
    % disp(Snew)
    F_cov_inv = inv(F2new);
    F = Fnew;

end

Final = F(:,1:L);
[LL L] = size(Final);  
Finalbar = sum(Final,2)/L
F3 = zeros(3);
for i=1:L
    F3 = F3 + (Final(:,i) - Finalbar)*(Final(:,i) - Fbar)';
end
Final_cov = F3/(L-1)
F_cov_inv = inv(F_cov); %For PC1, PC3 & PC3