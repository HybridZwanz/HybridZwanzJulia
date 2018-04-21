Q0=import_data("./QFull.csv");
Q1=import_data("./QApprox.csv");
Q2=import_data("QApprox2.csv");
P0=import_data("./PFull.csv");
P1=import_data("./PApprox.csv");
P2=import_data("PApprox2.csv");

[px0,x0]=histcounts(Q0,100,'Normalization','probability');
[px2,x2]=histcounts(Q2,x0,'Normalization','probability');
[px1,x1]=histcounts(Q1,x0,'Normalization','probability');
%[px05,x05]=histcounts(Q05,x0,'Normalization','probability');
[pm0,m0]=histcounts(P0,100,'Normalization','probability');
[pm2,m2]=histcounts(P2,m0,'Normalization','probability');
[pm1,m1]=histcounts(P1,m0,'Normalization','probability');
%[pm05,m05]=histcounts(P05,m0,'Normalization','probability');

figure()
hold off
x=(x0(2:end)+x0(1:end-1))/2;
plot(x,px0,'r')
hold on
plot(x,px1,'b')
plot(x,px2,'g')
%plot(x,px05,'g')

%dlmwrite("x.csv",[x;px0;px1;px2]')

figure()
hold off
m= (m0(2:end)+m0(1:end-1))/2;
plot(m,pm0,'r')
hold on
plot(m,pm1,'b')
plot(m,pm2,'g')
%plot(m,pm05,'g')

%dlmwrite("m.csv",[m;pm0;pm1;pm2]')


for i=1:2000
msd05(i)=var(Q05-circshift(Q05,i));
end
for i=1:2000
msd0(i)=var(Q0-circshift(Q0,i));
end
for i=1:2000
msd1(i)=var(Q1-circshift(Q1,i));
end
for i=1:2000
msd2(i)=var(Q2-circshift(Q2,i));
end
for i=1:2000
msd00(i)=var(Q00-circshift(Q00,i));
end
for i=1:2000-10
msd2down(i)=min(msd2(i:i+10));
msd0down(i)=min(msd0(i:i+10));
msd1down(i)=min(msd1(i:i+10));
msd00down(i)=min(msd00(i:i+10));
msd05down(i)=min(msd05(i:i+10));
end
figure
plot(msd0down)
hold on
plot(msd1down)
plot(msd2down)
plot(msd00down)
plot(msd05down)