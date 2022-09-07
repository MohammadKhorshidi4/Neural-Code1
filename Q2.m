%%
clc
clear
tau=8;
dt=1e-2;
vth=15;
RI=20;
V(1)=0;
for i = 1:9999
    dv = (-V(i)*dt + RI*dt)/tau;
    V(i+1)= V(i) + dv;
    if V(i+1) > vth
        V(i)=70;
        V(i+1)=0;
    end
end
t=linspace(0,100,10000);
figure('WindowState','maximized')
plot(t,V)
title('Leaky Integrate and Fire Model With $\tau$ = 8ms','Interpreter','latex')
xlabel('$Time (ms)$','Interpreter','Latex')
ylabel('$Voltage (mv)$','Interpreter','Latex')
grid on
grid minor
saveas(gcf,'Fig14.png')
%%
clc
clear
tau=1;
dt=1e-2;
vth=15;
RI=20;
V(1)=0;
for i = 1:9999
    dv = (-V(i)*dt + RI*dt)/tau;
    V(i+1)= V(i) + dv;
    if V(i+1) > vth
        V(i)=70;
        V(i+1)=0;
    end
end
t=linspace(0,100,10000);
figure('WindowState','maximized')
plot(t,V)
title('Leaky Integrate and Fire Model With $\tau$ = 1ms','Interpreter','latex')
xlabel('$Time (ms)$','Interpreter','Latex')
ylabel('$Voltage (mv)$','Interpreter','Latex')
grid on
grid minor
saveas(gcf,'Fig15.png')


%%
t=0:0.001:0.015;

x = t.*exp(-t/0.0015);
x = normalize(x,'range')*20;





sp1 = poissonSpikeGen(200, 1, 1);
sp1= double(sp1);

pt1 = conv(sp1,x,'same');
t1 = 0:1e-3:1 - 1e-3;
% pt1=pt1;
tau=0.001;
R=0.001;
vth = 15;




for i = 2:1e3
    c1 = 0;
    s=0:i-1;

     rp1 = exp(-s/tau);
     rp2 = pt1(i-s);
     c1 = trapz(s,rp1.*rp2);

    v(i) = (R/tau) * c1;
    if i>1
        if v(i) > vth
            v(i-1)=70;
            v(i)=0;
            pt1(i+1)=0;
        end
    end
end



figure('WindowState','maximized')

plot(t1,v)
xlabel('Time(s)','Interpreter','latex')
ylabel('mV','Interpreter','latex')
title('Voltage','Interpreter','latex')
grid on
grid minor
saveas(gcf,'Fig16.png')

figure('WindowState','maximized')
plot(t1,pt1)
title('Current','Interpreter','latex')
xlabel('Time(s)','Interpreter','latex')
ylabel('Amper','Interpreter','latex')
grid on
grid minor
saveas(gcf,'Fig17.png')
%%
clc
clear

tau1=linspace(0.1,10,10)*1e-2;

t=0:0.001:0.015;
x = t.*exp(-t/0.0015);
x = normalize(x,'range')*60;

t1=0:0.0001:0.015;
x1 = t1.*exp(-t1/0.0015);
x1 = normalize(x1,'range');
figure('WindowState','maximized')
plot(t1,x1)
xlabel('Time (s)','Interpreter','latex')
ylabel('Current','Interpreter','latex')
title('Current kernel','Interpreter','latex')
grid on
grid minor


R=1;

dt=1e-3;



for k=1:100
    for j=1:10
        tau = tau1(j);
        f = 0;
        sp1=zeros(1,1000);
        while f<200
            V = zeros(1,1000);
            v1 = zeros(1,1000);
            sp2 = poissonSpikeGen(200, 1, 1);
            sp2 = double(sp2);
            sp1 = sp1+sp2;
            pt1 = conv(sp1,x,'same');

            for i = 1:999
                dv = (-V(i)*dt + R*pt1(i)*dt)/tau;
                V(i+1)= V(i) + dv;
                if V(i+1) > k
                    V(i)=70;
                    V(i+1)=0;
                    pt1(i+1)=0;
                    v1(i)=1;
               end
            end
            f=sum(v1);
        end
        a1 = find(v1);
%         a3 = a1(1:x:end);
        a2 = diff(a1)*1e-3;
        cv1(k,j) = std(a2)/(mean(a2)+1e-10);
        m1{k,j}=v1;
        m2{k,j}=V;
    end
end
  
cv1=(floor(cv1*10)/10);


tau1=linspace(0.1,10,10)*1e-2;

for i =1:10
    xti{i}=num2str(tau1(i)*100);
end



figure('WindowState','maximized')
contour(cv1)
xticks(1:10)
xticklabels(xti)
xlabel('$\tau (ms)$ ','Interpreter','latex')
ylabel('$N_{th}$','Interpreter','latex')
title('CV Contour Plot','Interpreter','latex')
colorbar
saveas(gcf,'Fig19.png')


%%
clc
clear

tau1=linspace(0.1,10,10)*1e-2;

t=0:0.001:0.015;
x = t.*exp(-t/0.0015);
x = normalize(x,'range')*60;




R=1;

dt=1e-3;



for k=1:100
    for j=1:10
        tau = tau1(j);
        f = 0;
        sp1=zeros(1,1000);
        while f<200
            V = zeros(1,1000);
            v1 = zeros(1,1000);
            sp2 = poissonSpikeGen(200, 1, 1);
            sp2 = double(sp2);
            s1 = sum(sp2);
            s2 = ceil(0.1*s1);
            s3 = randperm(s1,s2);
            s4 = find(sp2);
            sp2(s4(s3)) = -1;
            sp1 = sp1+sp2;
            pt1 = conv(sp1,x,'same');

            for i = 1:999
                dv = (-V(i)*dt + R*pt1(i)*dt)/tau;
                V(i+1)= V(i) + dv;
                if V(i+1) > k
                    V(i)=70;
                    V(i+1)=0;
                    pt1(i+1)=0;
                    v1(i)=1;
               end
            end
            f=sum(v1);
        end
        a1 = find(v1);
%         a3 = a1(1:x:end);
        a2 = diff(a1)*1e-3;
        cv1(k,j) = std(a2)/(mean(a2)+1e-10);
        m1{k,j}=v1;
        m2{k,j}=V;
    end
end
  
cv1=(floor(cv1*10)/10);


tau1=linspace(0.1,10,10)*1e-2;

for i =1:10
    xti{i}=num2str(tau1(i)*100);
end



figure('WindowState','maximized')
contour(cv1)
xticks(1:10)
xticklabels(xti)
xlabel('$\tau (ms)$ ','Interpreter','latex')
ylabel('$N_{th}$','Interpreter','latex')
title('CV Contour Plot With 0.1 IPSP','Interpreter','latex')
colorbar
saveas(gcf,'Fig20.png')


%%
clc
clear

tau1=linspace(0.1,10,10)*1e-2;

t=0:0.001:0.015;
x = t.*exp(-t/0.0015);
x = normalize(x,'range')*60;




R=1;

dt=1e-3;



for k=1:100
    for j=1:10
        tau = tau1(j);
        f = 0;
        sp1=zeros(1,1000);
        while f<200
            V = zeros(1,1000);
            v1 = zeros(1,1000);
            sp2 = poissonSpikeGen(200, 1, 1);
            sp2 = double(sp2);
            s1 = sum(sp2);
            s2 = ceil(0.4*s1);
            s3 = randperm(s1,s2);
            s4 = find(sp2);
            sp2(s4(s3)) = -1;
            sp1 = sp1+sp2;
            pt1 = conv(sp1,x,'same');

            for i = 1:999
                dv = (-V(i)*dt + R*pt1(i)*dt)/tau;
                V(i+1)= V(i) + dv;
                if V(i+1) > k
                    V(i)=70;
                    V(i+1)=0;
                    pt1(i+1)=0;
                    v1(i)=1;
               end
            end
            f=sum(v1);
        end
        a1 = find(v1);
%         a3 = a1(1:x:end);
        a2 = diff(a1)*1e-3;
        cv1(k,j) = std(a2)/(mean(a2)+1e-10);
        m1{k,j}=v1;
        m2{k,j}=V;
    end
end
  
cv1=(floor(cv1*10)/10);


tau1=linspace(0.1,10,10)*1e-2;

for i =1:10
    xti{i}=num2str(tau1(i)*100);
end



figure('WindowState','maximized')
contour(cv1,'ShowText','on')
xticks(1:10)
xticklabels(xti)
xlabel('$\tau (ms)$ ','Interpreter','latex')
ylabel('$N_{th}$','Interpreter','latex')
title('$CV Contour Plot With 0.4 IPSP$','Interpreter','latex')
colorbar
saveas(gcf,'Fig21.png')
%%
clc
clear


R=1;
tau=.1;

t=0:0.001:0.015;
x = t.*exp(-t/0.0015);
x = normalize(x,'range')*10;

dt=1e-3;

ips = 1:49;
ips = ips./100;


for k = 1:length(ips)
    f = 0;
    sp1 = zeros(1,1000);
    while f<200
        V = zeros(1,1000);
        v1 = zeros(1,1000);
        sp2 = poissonSpikeGen(200, 1, 1);
        sp2 = double(sp2);
        s1 = sum(sp2);
        s2 = ceil((ips(k))*s1);
        s3 = randperm(s1,s2);
        s4 = find(sp2);
        sp2(s4(s3)) = -1;
        sp1 = sp1+sp2;
        pt1 = conv(sp1,x,'same');
        for i = 1:999
            dv = (-V(i)*dt + R*pt1(i)*dt)/tau;
            V(i+1)= V(i) + dv;
            if V(i+1) > 10
                V(i)=70;
                V(i+1)=0;
                pt1(i+1)=0;
                v1(i)=1;
            end
        end
        f=sum(v1);
    end
    a1 = find(v1);
%         a3 = a1(1:x:end);
    a2 = diff(a1)*1e-3;
    cv1(k) = std(a2)/(mean(a2)+1e-10);
    m1{k}=v1;
    m2{k}=V;
    pt1 = conv(sp1,x,'same');
    m3{k}=pt1;
end




figure('WindowState','maximized')
plot(ips.*100,cv1)
xlim([0,50])
xlabel('IPSP Percentage','Interpreter','latex')
ylabel('CV','Interpreter','latex')
title('Effects of IPSP Perentage Changes on CV','Interpreter','latex')
grid on
grid minor
saveas(gcf,'Fig22.png')



%%
clc
clear

wi = 0.0015:0.0001:0.0105;
mag = 10:90;





R=1;

dt=1e-3;

tau = 1;


for k=1:length(mag)
    for j=1:length(wi)

        t=0:0.001:0.015;
        
        tp = wi(j);

        x = t.*exp(-t/tp);
        
        x = normalize(x,'range') * mag(k);

        f = 0;
        sp1=zeros(1,1000);
        while f<200
            V = zeros(1,1000);
            v1 = zeros(1,1000);
            sp2 = poissonSpikeGen(200, 1, 1);
            sp2 = double(sp2);
%             s1 = sum(sp2);
%             s2 = ceil(0.4*s1);
%             s3 = randperm(s1,s2);
%             s4 = find(sp2);
%             sp2(s4(s3)) = -1;
            sp1 = sp1+sp2;
            pt1 = conv(sp1,x,'same');

            for i = 1:999
                dv = (-V(i)*dt + R*pt1(i)*dt)/tau;
                V(i+1)= V(i) + dv;
                if V(i+1) > 10
                    V(i)=70;
                    V(i+1)=0;
                    pt1(i+1)=0;
                    v1(i)=1;
               end
            end
            f=sum(v1);
        end
        a1 = find(v1);
%         a3 = a1(1:x:end);
        a2 = diff(a1)*1e-3;
        cv1(k,j) = std(a2)/(mean(a2)+1e-10);
        m1{k,j}=v1;
        m2{k,j}=V;
    end
end

figure('WindowState','maximized')
imagesc(cv1)
yticks(1:10:81)
xticks(1:10:91)

for i=1:9
    xti{i} = num2str(mag((i-1)*10+1));
end
for i=1:10
    yti{i} = num2str(wi((i-1)*10+1)*1e3);
end
xticklabels(yti)
yticklabels(xti)
xlabel('$\tau_{peak}(ms)$','Interpreter','latex')
ylabel('Current Kernel Magnitude','Interpreter','latex')
colorbar
title('CV changes over Currnet Kernel Magnitude and $\tau_{peak}$','Interpreter','latex')
saveas(gcf,'Fig23.png')



%%
clc
clear

wi = 0.015:0.001:0.055;
mag = 10:90;





R=1;

dt=1e-3;

tau = 2;


for k=1:length(mag)
    for j=1:length(wi)

        t=0:0.001:wi(j);
        
        tp = wi(j)/10;

        x = t.*exp(-t/tp);
        
        x = normalize(x,'range') * mag(k);

        f = 0;
        sp1=zeros(1,1000);
        while f<200
            V = zeros(1,1000);
            v1 = zeros(1,1000);
            sp2 = poissonSpikeGen(200, 1, 1);
            sp2 = double(sp2);
            
            sp1 = sp1+sp2;
            pt1 = conv(sp1,x,'same');

            for i = 1:999
                dv = (-V(i)*dt + R*pt1(i)*dt)/tau;
                V(i+1)= V(i) + dv;
                if V(i+1) > 10
                    V(i)=70;
                    V(i+1)=0;
                    pt1(i+1)=0;
                    v1(i)=1;
               end
            end
            f=sum(v1);
        end
        a1 = find(v1);
%         a3 = a1(1:x:end);
        a2 = diff(a1)*1e-3;
        cv1(k,j) = std(a2)/(mean(a2)+1e-10);
        m1{k,j}=v1;
        m2{k,j}=V;
    end
end

figure('WindowState','maximized')
imagesc(cv1)
yticks(1:10:81)
xticks(1:10:91)

for i=1:9
    xti{i} = num2str(mag((i-1)*10+1));
end
for i=1:5
    yti{i} = num2str(wi((i-1)*10+1)*1e3);
end
xticklabels(yti)
yticklabels(xti)
xlabel('Kernel Width (ms)','Interpreter','latex')
ylabel('Current Kernel Magnitude','Interpreter','latex')
colorbar
title('CV changes over Currnet Kernel Magnitude and Kernel Width','Interpreter','latex')
saveas(gcf,'Fig24.png')
%%
clc
clear




sp1 = zeros(1,1000);



for k = 1:100
    sp2 = poissonSpikeGen(50, 1, 1);
    sp2 = double(sp2);
    sp1 = sp1+sp2;
end

for i = 1:150
    for j=1:100
        v1 = zeros(1,floor(1000/i));
        for k=1:floor(1000/i)
            s1=sum(sp1(1,(k-1)*i+1:(k)*i));
            if s1 >= j
                v1(k) = 1;
            end
        end
        a1 = find(v1);
%         a3 = a1(1:x:end);
        a2 = diff(a1)*1e-3;
        cv1(i,j) = std(a2)/(mean(a2)+1e-10);
        m1{i,j}=v1;
%         m2{k,j}=V;
    end
end
        

figure('WindowState','maximized')
imagesc(cv1)

xticks(5:5:100)
c1=1;
for i =5:5:100
    xti{c1}=num2str(i/100);
    c1=c1+1;
end
xticklabels(xti)
ylabel('D (ms)','Interpreter','latex')
xlabel('$\frac{N}{M}$','Interpreter','latex')
title('CV Changes over D and $\frac{N}{M}$','Interpreter','latex')


colorbar
saveas(gcf,'Fig25.png')

figure('WindowState','maximized')
imagesc(cv1(1:20,:))

xticks(5:5:100)
c1=1;
for i =5:5:100
    xti{c1}=num2str(i/100);
    c1=c1+1;
end
xticklabels(xti)
ylabel('D (ms)','Interpreter','latex')
xlabel('$\frac{N}{M}$','Interpreter','latex')
title('CV Changes over D and $\frac{N}{M}$','Interpreter','latex')


colorbar
saveas(gcf,'Fig26.png')

%%

clc
clear




sp1 = zeros(1,1000);



for k = 1:100
    sp2 = poissonSpikeGen(50, 1, 1);
    sp2 = double(sp2);
    s1 = sum(sp2);
    s2 = ceil(0.4*s1);
    s3 = randperm(s1,s2);
    s4 = find(sp2);
    sp2(s4(s3)) = -1;
    sp1 = sp1+sp2;
end

for i = 1:150
    for j=1:100
        v1 = zeros(1,floor(1000/i));
        for k=1:floor(1000/i)
            s1=sum(sp1(1,(k-1)*i+1:(k)*i));
            if s1 >= j
                v1(k) = 1;
            end
        end
        a1 = find(v1);
%         a3 = a1(1:x:end);
        a2 = diff(a1)*1e-3;
        cv1(i,j) = std(a2)/(mean(a2)+1e-10);
        m1{i,j}=v1;
%         m2{k,j}=V;
    end
end
        

figure('WindowState','maximized')
imagesc(cv1)


ylabel('D (ms)','Interpreter','latex')
xlabel('$N_{net}$','Interpreter','latex')
title('CV Changes over D and $N_{net}$ with IPSP = 0.3','Interpreter','latex')
colorbar
saveas(gcf,'Fig27.png')


