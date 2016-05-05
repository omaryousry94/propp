clc
close all
clear all
rm=0.25;
cp=1005;
gamma=1.4;
R=287;
beta1dash=48.1;
beta2dash=12.9;
istar=2;
segma=0.94;
u=260;
m_dot=40;
pt0=1e5;
tt0=298;
beta1=istar+beta1dash;
theta=beta1dash-beta2dash;
beta2=((0.23*theta/sqrt(segma))+beta2dash)/(1-(theta/500/sqrt(segma)));
alpha2=beta1; alpha1=beta2;
cx1=(u/(tand(alpha1)+tand(beta1)));
c1=cx1/cosd(alpha1);
c0=cx1;
t0=tt0-(c0^2/2/cp);
m0=c0/sqrt(gamma*R*t0);
A0=((sqrt(gamma)*m0*(1+0.2*m0^2)^-3)*(pt0/m_dot/sqrt(R*tt0)))^-1;
conv_mfp=1/40;

%%
tt1=tt0;
y_igv=0.07;
t1=tt1-c1^2/(2*cp);
m1=c1/sqrt((gamma*R*t1));
syms A1sym
eq=(y_igv*(1-(1+0.2*m1^2)^-3.5)+1)^-1;
eq2=m_dot*sqrt(R*tt1)/(A1sym*eq*pt0*cosd(alpha1))-sqrt(gamma)*m1*(1+0.2*m1^2)^-3==0;
A1=double(solve(eq2,A1sym));
pt1=eq*pt0
%%
pt1_v=zeros(1,16);
tt1_v=zeros(1,16);
pist=zeros(1,8);
taust=zeros(1,8);
A=zeros(1,16);
t1_v=zeros(1,16);


pt1_v(1)=pt1;
tt1_v(1)=tt1;
pic=1;

for i=1:8
   
delta_c_theta=u-2*c1*sind(alpha1);
tt2=u*(delta_c_theta/cp)+tt1_v(i);

c2=cx1/cosd(beta1);
w1=c2;
w2=c1;
t2=tt2-c2^2/2/cp;
ttrel2=t2+w2^2/2/cp;
ttrel1=ttrel2;
ptrel1=pt1_v(i)*(ttrel1/tt1_v(i))^3.5;
t1_v(i)=tt1_v(i)-c1^2/(2*cp);
p1=pt1_v(i)*(t1_v(i)/tt1_v(i))^3.5;

w=0.07;
ptrel2=ptrel1-w*(ptrel1-p1);
pt2=ptrel2*(tt2/ttrel2)^3.5;
p2=ptrel2*(t2/ttrel2)^3.5;

m2=c2/sqrt(gamma*R*t2);

A(2*i-1)=((sqrt(gamma)*m2*(1+0.2*m2^2)^-3)*(pt2*cosd(alpha2)/m_dot/sqrt(R*tt2)))^-1;

tt3=tt2;
c3=c1;
alpha3=alpha1;
t3=tt3-c3^2/2/cp;
m3=c3/sqrt(gamma*R*t3);
pt3=pt2-w*(pt2-p2);
pt1_v(i+1)=pt3;
pist(i)=pt1_v(i+1)/pt1_v(i);
taust(i)=tt3/tt1;
tt1_v(i+1)=tt3;
A(2*i)=((sqrt(gamma)*m3*(1+0.2*m3^2)^-3)*(pt3*cosd(alpha3)/m_dot/sqrt(R*tt3)))^-1;

pic=pic*pist(i)
end
Atot=[A0 A1 A]
h=Atot/(2*pi*rm)

mfp0=(sqrt(R*tt3)/A(16)/pt3)*m_dot;

%% off design


epsiR=beta1-beta2;
epsiS=alpha1-alpha2;
deltaR=beta2-beta2dash;
alpha2dash=beta1dash;
alpha1dash=beta2dash;
deltaS=alpha1-alpha1dash;

Mrel2off=0.2; %initial Guess for Fzero Function Solution
M3off=0.2; %initial Guess for Fzero Function Solution
nm=200;
nn=8;
m_dot_v=linspace(0.3*m_dot,1.3*m_dot,nm);
N=linspace(6000,10500,nn)
%nc=8;
u_v=2*pi/60*N*rm;
m_dotrel_op=NaN(1,length(nn));
pi_op=NaN(1,length(nn));
%m_dot_v=40;
%u_v=260;
mat1=NaN(nn,nm);
mat2=NaN(nn,nm);
mat3=NaN(nn,nm);
inc=NaN(1,16);


% count=nc;
% no=0;
for t=1:length(u_v)
%for N
%for M_dot
for j=1:length(m_dot_v)
m_dotoff=m_dot_v(j);
tt1off=tt0; 
i=0;
Woff_igv=.56944*(i.^2)+.07;
%syms m1off
%eq=(Woff_igv*(1-(1+0.2*m1off^2)^-3.5)+1)^-1;
%eq2=m_dotoff*sqrt(R*tt1off)/(A1*eq*pt0*cosd(alpha1))-sqrt(gamma)*m1off*(1+0.2*m1off^2)^-3==0;
%M1off=double(solve(eq2,m1off));

M1off=fzero(@(m1off) (m_dotoff*sqrt(R*tt1off)/(A1*((Woff_igv*(1-(1+0.2*m1off^2)^-3.5)+1)^-1)*pt0*cosd(alpha1))-sqrt(gamma)*m1off*(1+0.2*m1off^2)^-3),0.3);
pt1off=(Woff_igv.*(1-(1+0.2.*M1off.^2).^-3.5)+1).^-1.*pt0;
p1off=pt1off./(1+0.2.*M1off.^2).^3.5;
tt1off=tt0;

alpha1off=alpha1;
uoff=u_v(t); %CHANGE AFTER FOR LOOP

tt1off_v=zeros(1,9);
pt1off_v=zeros(1,9);

alpha1off_v=zeros(1,9);

M1off_v=zeros(1,9);
pist_off=zeros(1,9);
taust_off=zeros(1,9);
tt1off_v(1)=tt1off;
pt1off_v(1)=pt1off;
alpha1off_v(1)=alpha1off;
M1off_v(1)=M1off;
picoff=1;
taucoff=1;
out=0;
%loop stage
for i=1:8
t1off=tt1off_v(i)./(1+0.2.*M1off_v(i).^2);
p1off=pt1off_v(i)./(1+0.2.*M1off_v(i).^2).^3.5;
c1off=M1off_v(i).*sqrt(gamma*R.*t1off);
cx1off=c1off.*cosd(alpha1off_v(i));
beta1off=atand(uoff./cx1off-tand(alpha1off_v(i)));
w1off=cx1off./cosd(beta1off);
ttrel1off=t1off+(w1off.^2./(2*cp));
ptrel1off=p1off.*(ttrel1off./t1off).^3.5;
ttrel2off=ttrel1off;
ir=beta1off-beta1dash;
k=(ir-istar)./epsiR;
Woff=.56944*(k.^2)+.07;
D=.9375*(k.^2);
deltaRoff=D.*epsiR+deltaR;
beta2off=deltaRoff+beta2dash;
ptrel2off=ptrel1off-Woff*(ptrel1off-p1off);
relmfp=m_dotoff*sqrt(R*ttrel2off)/(A(2*i-1)*ptrel2off*cosd(beta2off));
if relmfp>=0.684
    out=1;
    break
    
end
Mrel2off=fzero(@(mrel2off) (m_dotoff*sqrt(R*ttrel2off)/(A(2*i-1)*ptrel2off*cosd(beta2off))-sqrt(gamma)*mrel2off*(1+0.2*mrel2off^2)^-3),Mrel2off);

t2off=ttrel2off/(1+0.2*Mrel2off^2);
p2off=ptrel2off/(ttrel2off/t2off)^3.5;
w2off=Mrel2off*sqrt(gamma*R*t2off);
cx2off=w2off*cosd(beta2off);
alpha2off=atand((uoff/cx2off)-tand(beta2off));
c2off=cx2off/cosd(alpha2off);
m2off=c2off/sqrt(gamma*R*t2off);

tt2off=t2off+c2off^2/(2*cp);
pt2off=p2off*(tt2off/t2off)^3.5;
%% stator
is=alpha2off-alpha2dash;
if t==1
    inc(2*i-1)=ir;
    inc(2*i)=is;
    
end
ks=(is-istar)/epsiS;
Wsoff=.56944*(ks.^2)+.07;
if ks<0
    Ds=0;
else
Ds=.9375*(ks.^2);
end
deltasoff=Ds.*epsiS+deltaS;
pt3off=pt2off-Wsoff*(pt2off-p2off);
alpha3off=deltasoff+alpha1dash;
tt3off=tt2off;
mfp=(m_dotoff*sqrt(R*tt3off)/(A(2*i)*pt3off*cosd(alpha3off)));
if mfp>=0.684
    out=1;
    break
    
end
M3off=fzero(@(m3off) (m_dotoff*sqrt(R*tt3off)/(A(2*i)*pt3off*cosd(alpha3off))-sqrt(gamma)*m3off*(1+0.2*m3off^2)^-3),M3off);
t3off=tt3off/(1+0.2*M3off^2);
c3off=M3off*sqrt(gamma*R*t3off);
tt1off_v(i+1)=tt3off;
pt1off_v(i+1)=pt3off;
alpha1off_v(i+1)=alpha3off;
M1off_v(i+1)=M3off;
pist_off(i)=pt1off_v(i+1)/pt1off_v(i);
taust_off(i)=tt1off_v(i+1)/tt1off_v(i);
if pist_off(i)<1
    out=1;
    break
    
end

picoff=picoff*pist_off(i);
taucoff=taucoff*taust_off(i);
i

end
if out==1
    break
end

picoff_v(j)=picoff;
taucoff_v(j)=taucoff;
etacoff_v(j)=(picoff^(2/7)-1)/(taucoff-1);
%mat3(t,j)=etacoff_v(j);
j


if abs(mfp-mfp0)<=0.005;
m_dotrel_op(t)=m_dotoff*conv_mfp;
pi_op(t)=picoff_v(j);

end

end

figure(1)

[picoff_v_max(t),I] = max(picoff_v);
colormap('copper');
plot(m_dot_v(I:length(picoff_v))*conv_mfp,picoff_v(I:length(picoff_v)),'k','LineWidth', 1.47)
m_dot_v_max(t)=m_dot_v(I);
picoff_v_last(t)=picoff_v(j-1);
m_dot_v_last(t)=m_dot_v(j-1);

% if count==nc||t==nn;
% count=0;
% no=no+1;
% [picoff_v_max(no),I] = max(picoff_v);
% plot(m_dot_v(I:length(picoff_v))*conv_mfp,picoff_v(I:length(picoff_v)),'k')
% m_dot_v_max(no)=m_dot_v(I)
% picoff_v_last(no)=picoff_v(j-1);
% m_dot_v_last(no)=m_dot_v(j-1);
% 
% end
% count=count+1;



% mat3(t,j)=NaN;
% mat3(t,1:I_out)=NaN;
% mat1(t,I_out:j-1)=m_dot_v(I_out:j-1);
% mat2(t,I:j-1)=picoff_v(I:j-1);

xlim([0.2 1.2])
ylim([1 10])   
hold all
t
end
plot(m_dot_v_max*conv_mfp,picoff_v_max,'b','LineWidth', 1.6)
plot(m_dot_v_last*conv_mfp,picoff_v_last,'r','LineWidth', 1.6)
plot(m_dotrel_op,pi_op,'g','LineWidth', 1.6)


%% contour plotting
nm2=50;
for t=1:length(u_v)
%for N
%for M_dot
m_dot_limited=linspace(m_dot_v_max(t),m_dot_v_last(t),nm2);
for j=1:length(m_dot_limited)
m_dotoff=m_dot_limited(j);
tt1off=tt0; 
i=0;
Woff_igv=.56944*(i.^2)+.07;

M1off=fzero(@(m1off) (m_dotoff*sqrt(R*tt1off)/(A1*((Woff_igv*(1-(1+0.2*m1off^2)^-3.5)+1)^-1)*pt0*cosd(alpha1))-sqrt(gamma)*m1off*(1+0.2*m1off^2)^-3),0.3);
pt1off=(Woff_igv.*(1-(1+0.2.*M1off.^2).^-3.5)+1).^-1.*pt0;
p1off=pt1off./(1+0.2.*M1off.^2).^3.5;
tt1off=tt0;

alpha1off=alpha1;
uoff=u_v(t); %CHANGE AFTER FOR LOOP

tt1off_v=zeros(1,9);
pt1off_v=zeros(1,9);

alpha1off_v=zeros(1,9);

M1off_v=zeros(1,9);
pist_off=zeros(1,9);
taust_off=zeros(1,9);

tt1off_v(1)=tt1off;
pt1off_v(1)=pt1off;
alpha1off_v(1)=alpha1off;
M1off_v(1)=M1off;
picoff=1;
taucoff=1;
out=0;
%loop stage
for i=1:8
t1off=tt1off_v(i)./(1+0.2.*M1off_v(i).^2);
p1off=pt1off_v(i)./(1+0.2.*M1off_v(i).^2).^3.5;
c1off=M1off_v(i).*sqrt(gamma*R.*t1off);
cx1off=c1off.*cosd(alpha1off_v(i));
beta1off=atand(uoff./cx1off-tand(alpha1off_v(i)));
w1off=cx1off./cosd(beta1off);
ttrel1off=t1off+(w1off.^2./(2*cp));
ptrel1off=p1off.*(ttrel1off./t1off).^3.5;
ttrel2off=ttrel1off;
ir=beta1off-beta1dash;

k=(ir-istar)./epsiR;
Woff=.56944*(k.^2)+.07;
D=.9375*(k.^2);
deltaRoff=D.*epsiR+deltaR;
beta2off=deltaRoff+beta2dash;
ptrel2off=ptrel1off-Woff*(ptrel1off-p1off);
relmfp=m_dotoff*sqrt(R*ttrel2off)/(A(2*i-1)*ptrel2off*cosd(beta2off));
% if relmfp>=0.684
%     out=1;
%     break
%     
% end
Mrel2off=fzero(@(mrel2off) (m_dotoff*sqrt(R*ttrel2off)/(A(2*i-1)*ptrel2off*cosd(beta2off))-sqrt(gamma)*mrel2off*(1+0.2*mrel2off^2)^-3),0.3);

t2off=ttrel2off/(1+0.2*Mrel2off^2);
p2off=ptrel2off/(ttrel2off/t2off)^3.5;
w2off=Mrel2off*sqrt(gamma*R*t2off);
cx2off=w2off*cosd(beta2off);
alpha2off=atand((uoff/cx2off)-tand(beta2off));
c2off=cx2off/cosd(alpha2off);
m2off=c2off/sqrt(gamma*R*t2off);

tt2off=t2off+c2off^2/(2*cp);
pt2off=p2off*(tt2off/t2off)^3.5;
%% stator
is=alpha2off-alpha2dash;
ks=(is-istar)/epsiS;
Wsoff=.56944*(ks.^2)+.07;
if ks<0
    Ds=0;
else
Ds=.9375*(ks.^2);
end
deltasoff=Ds.*epsiS+deltaS;
pt3off=pt2off-Wsoff*(pt2off-p2off);
alpha3off=deltasoff+alpha1dash;
tt3off=tt2off;
mfp=(m_dotoff*sqrt(R*tt3off)/(A(2*i)*pt3off*cosd(alpha3off)));
% if mfp>=0.684
%     out=1;
%     break
%     
% end
M3off=fzero(@(m3off) (m_dotoff*sqrt(R*tt3off)/(A(2*i)*pt3off*cosd(alpha3off))-sqrt(gamma)*m3off*(1+0.2*m3off^2)^-3),0.3);
t3off=tt3off/(1+0.2*M3off^2);
c3off=M3off*sqrt(gamma*R*t3off);
tt1off_v(i+1)=tt3off;
pt1off_v(i+1)=pt3off;
alpha1off_v(i+1)=alpha3off;
M1off_v(i+1)=M3off;
pist_off(i)=pt1off_v(i+1)/pt1off_v(i);
taust_off(i)=tt1off_v(i+1)/tt1off_v(i);


picoff=picoff*pist_off(i);
taucoff=taucoff*taust_off(i);
i

end


picoff_v(j)=picoff;
taucoff_v(j)=taucoff;
etacoff_v(j)=(picoff^(2/7)-1)/(taucoff-1);
mat3(t,j)=etacoff_v(j);
j


end

 mat1(t,1:nm2)=m_dot_limited(1:nm2);
 mat2(t,1:nm2)=picoff_v(1:nm2);

t
end
[C1,h1]=contour(mat1*conv_mfp,mat2,mat3,[0.78:0.01:0.85],'LineWidth', 1.3)
%clabel(C1,h1)
set(h1,'showText','on','TextStep',get(h1,'LevelStep')*3)
colormap('winter');
figure(2)
surge_margin=(picoff_v_max-pi_op)./pi_op;
plot(m_dotrel_op,surge_margin)
grid on
xlabel('MFP_r_e_l')
ylabel('Surge Margin')
figure(3)
plot([1:1:16],inc)