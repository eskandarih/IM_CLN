tic
clear
format long e

files = fopen('./s85m17h.dat','r');
filem = fopen('./m85m17h.dat','r');
lambda = fscanf(files,'%e');

m=lambda(1)-3;  %2*harmonic
NN=4;           %段数
s=0.5;          %滑り
freq=50;        %周波数
w=2*pi*freq;    %角周波数  
Tm=zeros(2,2);  %変換局所行列
T=sparse(m,m);  %変換行列
per=3;                 %周期数
nn=360;                %時間分割数
nmax=nn+1;             %時間分割+1
tmax=per/freq;         %最大時間
dt=tmax/nn;            %時間幅
r=15.3e-3;
pp=2;                   %極対数

filev = fopen('v050s05n4p85m17h.dat','w');
filep = fopen('p050s05n4p85m17h.dat','w');

%%stator
Ls0=zeros(m+3,m+3);
for j=1:1:m+3 
  for k=1:1:m+3 
      Ls(j,k)=lambda((j-1)*(m+3)+k+2);
  end
end
L00=Ls(1:3,1:3);
L01=Ls(1:3,4:m+3);
L11=Ls(4:m+3,4:m+3);
L10=L01';
fclose(files);

%%mover
Lra=zeros(m*NN,m*NN);
Rra=zeros(m*NN,m*NN);
Dr =eye(m*NN);
Lr=zeros(m,m,NN);
Rr=zeros(m,m,NN);

lambda = fscanf(filem,'%e');
fclose(filem);
ll=2;
for i=1:1:NN
  for j=1:1:m 
    for k=1:1:m 
      Lr(j,k,i)=lambda((j-1)*m+k+ll);
      Rr(j,k,i)=lambda((j-1)*m+k+ll+m*m+2);
    end
  end
  ll=ll+2*m*m+4;
end
for i=1:1:NN
  Rr(:,:,i)=inv(Rr(:,:,i));
  Lra((i-1)*m+1:i*m,(i-1)*m+1:i*m)=Lr(:,:,i);
  Rra((i-1)*m+1:i*m,(i-1)*m+1:i*m)=Rr(:,:,i);
end

for i=1:1:NN-1
  Dr((i-1)*m+1:i*m,i*m+1:(i+1)*m)=-eye(m);
end
Drinv=inv(Dr);

Lra=Dr*Lra;
Rra=Rra*Drinv';

t=0;
ii=zeros(m*NN+m,1);                %電流初期値I'2n-1
In=zeros(m*NN+m,nmax);             %I'2n-1,I
I00=zeros(3,1); 
I0 =zeros(3,1); 

  I00=[cos(0.0); cos(-2*pi/3); cos(-4*pi/3)];   %一次電流の初期値 
  ii(m*NN+1:m*NN+m)=-inv(L11)*L10*I00;          %電流初期値I
  In(:,1)=ii;

%時間依存微分方程式を解く
%t→t+dt
for j=1:1:nn 
%for j=1:1:nmax 
  
  L=sparse(m*NN+m,m*NN+m);
  R=sparse(m*NN+m,m*NN+m);
  A=zeros(m*NN+m,m*NN+m);
  B=zeros(m*NN+m,1);
  
  I0=[cos(w*(t+dt)); cos(w*(t+dt)-2*pi/3); cos(w*(t+dt)-4*pi/3)];   %一次電流 
%  I0=[cos(w*t); cos(w*t-2*pi/3); cos(w*t-4*pi/3)];   %一次電流 
  
  %変換行列T
  for k=1:2:m
    Tm(:,:)=[cos(k*(1-s)*w*(t+dt)) sin(k*(1-s)*w*(t+dt));
            -sin(k*(1-s)*w*(t+dt)) cos(k*(1-s)*w*(t+dt))];
%    Tm(:,:)=[cos(k*(1-s)*w*t) sin(k*(1-s)*w*t);
%            -sin(k*(1-s)*w*t) cos(k*(1-s)*w*t)];
    T(k:k+1,k:k+1)=Tm;
  end

  %インダクタの合成行列L
  L(1:m*NN,1:m*NN)=Lra;
  L(m*NN+1:m*NN+m,1:m)=inv(T)*Lr(:,:,1);
  L(m*NN+1:m*NN+m,m*NN+1:m*NN+m)=-L11;
  
  %抵抗の合成行列R
  R(1:m*NN,1:m*NN)=Rra;
  for i=1:1:NN
    R((i-1)*m+1:i*m,m*NN+1:m*NN+m)=Rr(:,:,i)*T;
  end
  
  %入力電流作用ベクトルB
  B=L*ii;
  B(m*NN+1:m*NN+m)=L10*I0;
  
  A=L+dt*R;
  
  %電流I'2n-1,I
  ii=inv(A)*B;
  t=t+dt;
  In(:,j+1)=ii;
%  In(:,j)=ii;

end

In=In';
t=linspace(0,tmax,nmax);
%disp(In);

I  =zeros(m,nmax);
Ph =zeros(m,nmax);
Ph0=zeros(3,nmax);
V0 =zeros(3,nmax);
V  =zeros(m,nmax);
P0 =zeros(1,nmax);
P  =zeros(1,nmax);
TT =zeros(m,nmax);
Tx =zeros(1,nmax);
Ph00=zeros(3,1);
Pho=zeros(m,1);
Io =zeros(m,1);

t=0;

%I',I,Ph0,Phの算出
for n=1:1:nmax
  I0=[cos(w*t); cos(w*t-2*pi/3);  cos(w*t-4*pi/3)];   %一次電流
  II=In(n,m*NN+1:m*NN+m)';  %電流I
  I(:,n)  =II; 
  Ph0(:,n)=L00*I0+L01*II;               %一次磁束PHI0
  Ph(:,n) =L10*I0+L11*II;               %空隙磁束PHI   
  V0(:,n) =(Ph0(:,n)-Ph00(:))/dt;       %一次電圧V0
  V(:,n)  =(Ph(:,n)-Pho(:))/dt;         %空隙電圧V
  P0(1,n) =V0(:,n)'*(I0+I00)/2;         %一次入力
  P(1,n)  =-V(:,n)'*(I(:,n)+Io)/2;      %二次入力(空隙部電力)
  
  t=t+dt;
  Ph00(:)=Ph0(:,n);
  Pho(:) =Ph(:,n);
  I00    =I0;
  Io(:)  =I(:,n);
end
%disp(Ph0');

t=dt;          %n=2からなのでtはdtから


%トルクの算出
%最終周期の平均を取る

for k=1:2:m
  for n=1:1:nmax
    TT(k,n)=-Ph(k,n)*I(k+1,n)+Ph(k+1,n)*I(k,n);
  end
end  

for k=1:2:m
  Tx(1,:)=Tx(1,:)+TT(k,:)*k;
end
Tx=pp*Tx;
for n=1:1:nmax
  v(1)=(n-1)*dt;
  v(2:4)= Ph0(:,n);
  v(5:7)= V0(:,n);
  fprintf(filev,'%e %e %e %e %e %e %e\n',v);
  u(1)=(n-1)*dt;
  u(2)=Tx(1,n);
  u(3)=P0(1,n);
  u(4)=P(1,n);
  u(5)=u(2)*(1-s)*w;
  fprintf(filep,'%e %e %e %e %e\n',u);
end

%disp(I(:,nn*(per-1)/per+1:nn)');
%disp(Ph(:,nn*(per-1)/per+1:nn)');

disp(Tx');
tx=mean(Tx(:,nn*(per-1)/per+1:nn));
disp(tx'); 
fclose(filev);
fclose(filep);

