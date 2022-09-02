%%
clc
clear
[spikeMat, tVec] = poissonSpikeGen(100, 1, 100);
figure('WindowState','maximized')
plotRaster(spikeMat, tVec)
xlabel('Time(S)','Interpreter','latex')
ylabel('Number of Trials','Interpreter','latex')
saveas(gcf,'Fig1.png')
%%
clc
clear
for i=1:1000
    [spikeMat, tVec] = poissonSpikeGen(100, 1, 1);
    sp=double(spikeMat);
    sp1(i)=sum(sp,["all"]);
end
figure('WindowState','maximized')
h = histogram(sp1,'Normalization','pdf');
% h.values
hold on

x=1:200;
y=poisspdf(x,100);
plot(y,'Color','red')
xlabel('rate','Interpreter','latex')
ylabel('Probabilty','Interpreter','latex')
legend('Spikes count pdf','Poisson Distrobution with \lambda = 100')
saveas(gcf,'Fig2.png')


%%
clc
clear
% t=0:0.0001:50;
% ex1=exp(-1000*t);
[spikeMat, tVec] = poissonSpikeGen(100, 1, 1);
sp1=find(spikeMat);
sp2=diff(sp1*0.001);
figure('WindowState','maximized')
h = histogram(sp2,'Normalization','probability','NumBins',5);
xlim([0,0.1])
hold on
x1 = 0.05:0.001:.2;
y = 100*exp(-100*x1);
plot(x1-0.05,y,'Color','red')
xlabel('$\Delta$t','Interpreter','latex')
ylabel('Probability','Interpreter','latex')
legend('ISI pdf','Exponential Distrobution pdf with \lambda=100')

%%

clc
clear
[spikeMat, tVec] = poissonSpikeGen(100, 1, 1000);
c1 = 0;
c2 = 0;
c3 = 0;
c4 = 0;
sp1 = zeros(1000,1000);
sp2 = zeros(1000,1000);
sp3 = zeros(1000,1000);
sp4 = zeros(1000,1000);

for i = 1:1000
    for j = 1:1000
        if spikeMat(i,j)
            c1 = c1+1;
            c2 = c2+1;
            c3 = c3+1;
            c4 = c4+1;
            if c1 == 2
                sp1(i,j) = 1;
                c1 = 0 ;
            end
            if c2 == 4
                sp2(i,j) = 1;
                c2 = 0;
            end
            if c3 == 6
                sp3(i,j) = 1;
                c3 = 0;
            end
            if c4 == 10
                sp4(i,j) = 1;
                c4 = 0;
            end
        end
    end
end

s1 = sum(sp1');
s2 = sum(sp2');
s3 = sum(sp3');
s4 = sum(sp4');

figure('WindowState','maximized')
h = histogram(s1,'Normalization','pdf');
hold on

x=1:200;
y=poisspdf(x,100);
plot(y,'Color','red')

xticks(0:5:200)
xlabel('rate','Interpreter','latex')
ylabel('Probabilty','Interpreter','latex')
legend('Spikes count pdf for k = 2','Poisson Distrobution with \lambda = 100 with no integration')
saveas(gcf,'Fig4.png')

figure('WindowState','maximized')
h = histogram(s2,'Normalization','pdf');
hold on

x=1:200;
y=poisspdf(x,100);
plot(y,'Color','red')
xlabel('rate','Interpreter','latex')
ylabel('Probabilty','Interpreter','latex')
legend('Spikes count pdf for k = 4','Poisson Distrobution with \lambda = 100 with no integration')
saveas(gcf,'Fig5.png')
xticks(0:5:200)


figure('WindowState','maximized')
h = histogram(s3,'Normalization','pdf');
hold on

x=1:200;
y=poisspdf(x,100);
plot(y,'Color','red')
xlabel('rate','Interpreter','latex')
ylabel('Probabilty','Interpreter','latex')
legend('Spikes count pdf for k = 6','Poisson Distrobution with \lambda = 100 with no integration')
saveas(gcf,'Fig6.png')
xticks(0:5:200)


figure('WindowState','maximized')
h = histogram(s4,'Normalization','pdf');
hold on

x=1:200;
y=poisspdf(x,100);
plot(y,'Color','red')
xlabel('rate','Interpreter','latex')
ylabel('Probabilty','Interpreter','latex')
legend('Spikes count pdf for k = 10','Poisson Distrobution with \lambda = 100 with no integration')
saveas(gcf,'Fig7.png')
xticks(0:5:200)



%%
clc
clear
[spikeMat, tVec] = poissonSpikeGen(100, 1, 1000);
c1 = 0;
c2 = 0;
c3 = 0;
c4 = 0;
sp1 = zeros(1000,1000);
sp2 = zeros(1000,1000);
sp3 = zeros(1000,1000);
sp4 = zeros(1000,1000);

for i = 1:1000
    for j = 1:1000
        if spikeMat(i,j)
            c1 = c1+1;
            c2 = c2+1;
            c3 = c3+1;
            c4 = c4+1;
            if c1 == 2
                sp1(i,j) = 1;
                c1 = 0 ;
            end
            if c2 == 4
                sp2(i,j) = 1;
                c2 = 0;
            end
            if c3 == 6
                sp3(i,j) = 1;
                c3 = 0;
            end
            if c4 == 10
                sp4(i,j) = 1;
                c4 = 0;
            end
        end
    end
end
for i = 1:1000
    a = sp1(i,:);
    a1 = find(a);
    a2 = diff(a1); 

    a3 = sp2(i,:);
    a4 = find(a3);
    a5 = diff(a4);

    a6 = sp3(i,:);
    a7 = find(a6);
    a8 = diff(a7);

    a9 = sp4(i,:);
    a10 = find(a9);
    a11 = diff(a10);

    if i==1
        x1 = a2;
        x2 = a5;
        x3 = a8;
        x4 = a11;
    else
        x1 = [x1,a2];
        x2 = [x2,a5];
        x3 = [x3,a8];
        x4 = [x4,a11];
    end
end


figure('WindowState','maximized')
histogram(x1*0.001,'Normalization','probability')
hold on
erlang_fun = @(x, mu, k) x.^(k-1).*exp(-x/mu)/(mu.^k.*factorial(k-1));
x = 0:0.001:1;
x5 = 0.05:0.001:.2;
y = 10*exp(-100*x5);
plot(x5-0.05,y,'Color','red')
hold on
plot(x,erlang_fun(x,0.01,2)/1000,'Color','green','LineWidth',3)
xlabel('$\Delta$t','Interpreter','latex')
ylabel('Probability','Interpreter','latex')
legend('ISI pdf','Exponential Distrobution pdf with \lambda=100','Erlang Distrobution with k = 2 and \lambda=100')
saveas(gcf,'Fig8.png')


figure('WindowState','maximized')
histogram(x2*0.001,'Normalization','probability')
hold on
erlang_fun = @(x, mu, k) x.^(k-1).*exp(-x/mu)/(mu.^k.*factorial(k-1));
x = 0:0.001:1;
x5 = 0.05:0.001:.2;
y = 10*exp(-100*x5);
plot(x5-0.05,y,'Color','red')
hold on
plot(x,erlang_fun(x,0.01,4)/500,'Color','green','LineWidth',3)
xlabel('$\Delta$t','Interpreter','latex')
ylabel('Probability','Interpreter','latex')
legend('ISI pdf','Exponential Distrobution pdf with \lambda=100','Erlang Distrobution with k = 4 and \lambda=100')
saveas(gcf,'Fig9.png')

figure('WindowState','maximized')
histogram(x3*0.001,'Normalization','probability')
hold on
erlang_fun = @(x, mu, k) x.^(k-1).*exp(-x/mu)/(mu.^k.*factorial(k-1));
x = 0:0.001:1;
x5 = 0.05:0.001:.2;
y = 10*exp(-100*x5);
plot(x5-0.05,y,'Color','red')
hold on
plot(x,erlang_fun(x,0.01,6)/333,'Color','green','LineWidth',3)
xlabel('$\Delta$t','Interpreter','latex')
ylabel('Probability','Interpreter','latex')
legend('ISI pdf','Exponential Distrobution pdf with \lambda=100','Erlang Distrobution with k = 6 and \lambda=100')
saveas(gcf,'Fig10.png')

figure('WindowState','maximized')
histogram(x4*0.001,'Normalization','probability')
hold on
erlang_fun = @(x, mu, k) x.^(k-1).*exp(-x/mu)/(mu.^k.*factorial(k-1));
x = 0:0.001:1;
x5 = 0.05:0.001:.2;
y = 10*exp(-100*x5);
plot(x5-0.05,y,'Color','red')
hold on
plot(x,erlang_fun(x,0.01,10)/200,'Color','green','LineWidth',3)
xlabel('$\Delta$t','Interpreter','latex')
ylabel('Probability','Interpreter','latex')
legend('ISI pdf','Exponential Distrobution pdf with \lambda=100','Erlang Distrobution with k = 10 and \lambda=100')
saveas(gcf,'Fig11.png')

%% 
clc
clear

c1 = 0;
c2 = 0;
c3 = 0;
c4 = 0;

for i = 1:1000
    [spikeMat, tVec] = poissonSpikeGen(100, 1, 1);
    sp1 = double(spikeMat);
    
    a1 = find(sp1);
    a2 = diff(a1);

    a4 = a1(2:2:end);
    
    
    a5 = diff(a4);
    
    a7 = a1(4:4:end);
    
    a8 = diff(a7);
    
    a10 = a1(6:6:end);

    a11 = diff(a10);

    a13 = a1(10:10:end);

    a14 = diff(a13);

    cv1(i) = sqrt(var(a2))/mean(a2);
    cv2(i) = sqrt(var(a5))/mean(a5);
    cv3(i) = sqrt(var(a8))/mean(a8);
    cv4(i) = sqrt(var(a11))/mean(a11);
    cv5(i) = sqrt(var(a14))/mean(a14);
end
cv6=[cv1,cv2,cv3,cv4,cv5];

figure('WindowState','maximized')
for i = 1:5
    x = ((i-1)*1000 + 1):i*1000;
    bar(x,cv6(x))
    hold on
end
xticks(500:1000:4500);
xticklabels({'k=1';'k=2';'k=4';'k=6';'k=10'})
title('Coefficient of Variation','Interpreter','latex')
k1 = [1,2,4,6,10];
k2 = 1./sqrt(k1);
for i = 1:5
    x = ((i-1)*1000 + 1):i*1000;
    k3 = ones(1000,1);
    plot(x,k3*k2(i),'k--',LineWidth=1.5)
    text(250+(i-1)*1000,k2(i)+0.025,['$\downarrow$ Theoretical $C_v$ for k = ',num2str(k1(i))],...
        'Interpreter','latex')
    hold on
end
legend('Generated Trials C_v for k = 1','Generated Trials C_v for k = 4', ...
    'Generated Trials C_v for k = 6', ...
    'Generated Trials C_v for k = 10')
saveas(gcf,'Fig12.png')




%%
clc
clear
c3=0;
for i = 1000:-0.25:30
    c3=c3+1;
    sp101=0;
    spr77=0;
    [spikeMat, tVec] = poissonSpikeGen(i, 5, 1);
    sp1=find(spikeMat);
    for j=2:length(sp1)
        if sp1(j)-sp1(j-1)==1 && sp1(j-1) ~= 0
            sp1(j)=0;
        end
    end
    ind1=find(sp1~=0);
    sp11=sp1(ind1);
    sp1d=diff(sp11);
    cv1(c3)=sqrt(var(sp1d))/mean(sp1d);
    sp1=find(spikeMat);
    c1=0;
    for j=1:length(sp1)-1
        if sp1(j)~=0
            c1=c1+1;
        end
        if c1==4
            c1=0;
            if sp1(j+1)-sp1(j)==1
                sp1(j+1)=0;
            end
        end
    end 
    ind1=find(sp1~=0);
    sp12=sp1(ind1);
    sp44=sp12(4:4:end);
 
    sp1=find(spikeMat);
    c1=0;
    for j=1:length(sp1)-1
        if sp1(j)~=0
            c1=c1+1;
        end
        if c1==20
            c1=0;
            if sp1(j+1)-sp1(j)==1
                sp1(j+1)=0;
            end
        end
    end 
    ind1=find(sp1~=0);
    sp12=sp1(ind1);
    sp201=sp12(20:20:end);
    sp1=find(spikeMat);
    c1=0;
    c2=1;
    for j=1:length(sp1)-1
        if sp1(j)~=0
            c1=c1+1;
        end
        if c1==10
            c1=0;
            if sp1(j+1)-sp1(j)==1
                sp1(j+1)=0;
            end
        end
    end 
    ind1=find(sp1~=0);
    sp12=sp1(ind1);
    sp101=sp12(10:10:end);
    sp44d=diff(sp44);
    sp101d=diff(sp101);
    sp201d=diff(sp201);
    cv4(c3)=sqrt(var(sp44d))/mean(sp44d);
    cv10(c3)=sqrt(var(sp101d))/mean(sp101d);
    cv20(c3)=sqrt(var(sp201d))/mean(sp201d);
end
f12=1000./(1000:-0.25:30);

figure('WindowState','maximized')
plot(f12,cv1,'.')
hold on
yline(1)
xt = f12(10);
yt = 1.05;
str = 't_r_e_f = 0 & k=1';
text(xt,yt,str)
hold on
plot(f12,cv4,'.')
hold on
yline(0.5)
xt = f12(10);
yt = 0.55;
str = 't_r_e_f = 0 & k=4';
text(xt,yt,str)
hold on
plot(f12,cv10,'.')
hold on
yline(1/sqrt(10))
xt = f12(10);
yt = 1/sqrt(10)+.05;
str = 't_r_e_f = 0 & k=10';
text(xt,yt,str)
hold on
plot(f12,cv20,'.')
xlabel('$\bar{\Delta t} (ms)$','Interpreter','Latex')
ylabel('$C_v$','Interpreter','Latex')
hold on
yline(1/sqrt(20))
xt = f12(10);
yt = 1/sqrt(20)+.05;
str = 't_r_e_f = 0 & k=20';
text(xt,yt,str)
legend('t_r_e_f = 1ms & k=1','','t_r_e_f = 1ms & k=4','','t_r_e_f = 1ms & k=10',...
    '','t_r_e_f = 1ms & k=20')
saveas(gcf,'Fig13.png')