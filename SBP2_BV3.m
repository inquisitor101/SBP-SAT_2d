%m=11; %problemstorlek
%h=1/(m-1);

% Second order accurate SBP operators. Specify nummber of gripdpoints and h
% Note that h depend on the size of the computational domain (it is not
% always 1

%D2=HI(-A+BD)

e_1=zeros(m,1);e_1(1)=1;
e_m=zeros(m,1);e_m(m)=1;

H=(eye(m,m));H(1,1)=0.5;H(m,m)=0.5;
H=h*H;
HI=inv(H);

D1=((.5*diag(ones(m-1,1),1)-.5*diag(ones(m-1,1),-1)));
D1(1,1)=-1;D1(1,2)=1;D1(m,m-1)=-1;D1(m,m)=1;
D1(m,m-1)=-1;D1(m,m)=1;
D1=D1/h;

Q=H*D1 + 1/2*e_1*e_1' - 1/2*e_m*e_m';

D2=((diag(ones(m-1,1),1)+diag(ones(m-1,1),-1)-2*diag(ones(m,1),0)));
D2(1,1)=1;D2(1,2)=-2;D2(1,3)=1;
D2(m,m-2)=1;D2(m,m-1)=-2;D2(m,m)=1;
D2=D2/h^2;

S_U=[-3/2, 2, -1/2]/h;
S_1=zeros(1,m);
S_1(1:3)=S_U;
S_m=zeros(1,m);
S_m(m-2:m)=fliplr(-S_U);

M=-H*D2-e_1*S_1+e_m*S_m;
