%This Hamiltonian is designed for ultrathin TI films
%Spintexture fixed by LE QUY DUONG
%Written by LE QUY DUONG
%

function [ HH ] = genH(kx,ky,NLayer,Mz )
%%Based on Spin-orbit coupled 2D electron gas 
%%kx : GK
%%ky : GM
tpara=[
    -0.24
    12
    -18
    2.2
    0.39
    0.0
    -0.042
];
d=tpara(1);
m1=1/tpara(2);
m2=1/tpara(3);
v=tpara(4);
tz=tpara(5);
tz2=tpara(7);
lastQL=4*(NLayer-1)+1:4*NLayer;
firstQL=1:4;
mfield=[0 0 Mz];

kp=kx+1i*ky;
km=kx-1i*ky;
k2=kx^2+ky^2;

sigx=[0 0 1 0;0 0 0 1;1 0 0 0;0 1 0 0];
sigy=[0 0 -1i 0;0 0 0 -1i;1i 0 0 0;0 1i 0 0];
sigz=[1 0 0 0; 0 1 0 0; 0 0 -1 0; 0 0 0 -1];

T=[tz2 0 0 0;
   tz tz2 0 0;
   0 0 tz2 0;
   0 0 tz  tz2];

 Hp=[ k2/m1  d+k2/m2  v*kp    0;
  d+k2/m2  k2/m1     0    -v*kp;
  v*kp'  0  k2/m1    d+k2/m2;
  0   -v*kp'    d+k2/m2   k2/m1];

HH=kron(diag(linspace(1,1,NLayer),0),Hp) + kron(diag(linspace(1,1,NLayer-1),1),T) + kron(diag(linspace(1,1,NLayer-1),-1),T');
HH(lastQL,lastQL)=HH(lastQL,lastQL)+mfield(1)*sigx+mfield(2)*sigy+mfield(3)*sigz;

return