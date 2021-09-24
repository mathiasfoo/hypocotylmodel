function dC = proposed_RBLight_ODEmod_com_full(t,C)


%% Modified DeCaluwe's Model (Add Red and Blue Light and LHY inhibit PRR9)

global Ired Iblue

v1=4.8318;q1a=1.4266;q3a=8.9432;eta1=0.03;q4a=5.9277;eta2=0.0215;K1=0.1943;K2=1.6138;k1L=0.2866;k1D=0.213;
p1=0.8672;p1L=0.2378;d1=0.7843;q1b=3.575;q3b=5.5899;q4b=8.954;v2=1.6822;K3=2.2275;K4=0.40;K5=0.37;k2=0.35;
p2=0.7858;d2D=0.3712;d2L=0.2917;v3=1.113;K6=0.4944;K7=2.4087;k3=0.5819;
p3=0.6142;d3D=0.5026;d3L=0.5431;v4=2.5012;K8=0.3262;K9=1.7974;K10=1.1889;k4=0.925;
p4=1.126;de1=0.0022;de2=0.4741;de3=0.3765;de4=0.398;de5=0.0003;
Ap3=0.3868;Am7=0.5503;Ak7=1.125;q2=0.5767;kmpac=137;kd=7;
v5=0.1129;K11=0.3322;k5=0.1591;p5=0.5293;d5D=0.4404;d5L=5.0712;g1=0.001;g2=0.18;
K12=0.86;Bp4=0.4147;Bm8=0.7728;Bk8=0.1732;kmpbc=7162;Cp5=0.4567;Cm9=0.867;Ck9=0.3237;kmcc=13406;


if Ired ~= 0 || Iblue~= 0
    ThetaPhyA = 1;
else
    ThetaPhyA = 0;
end

if Ired ~= 0
    ThetaPhyB = 1;
else
    ThetaPhyB = 0;
end

if Iblue ~= 0
    ThetaCry1 = 1;
else
    ThetaCry1 = 0;
end

dC = zeros(18,1);

%LHY mRNA
dC(1) = (v1 + (q1a*(C(9)+C(16))*ThetaPhyA+q3a*(C(13)+C(17))*log(eta1*Ired + 1)*ThetaPhyB+q4a*(C(14)+C(18))*log(eta2*Iblue + 1)*ThetaCry1))/(1 + (C(4)/(K1))^2 + (C(6)/(K2))^2) - (k1L*ThetaPhyA + k1D*(1-ThetaPhyA))*C(1);

% LHY protein
dC(2) = (p1 + p1L*ThetaPhyA)*C(1) - d1*C(2);

% P97 mRNA
dC(3) = ((q1b*(C(9)+C(16))*ThetaPhyA+q3b*(C(13)+C(17))*log(eta1*Ired + 1)*ThetaPhyB+q4b*(C(14)+C(18))*log(eta2*Iblue + 1)*ThetaCry1)+ v2)*(1/(1 + (C(2)/(K3))^2 + (C(6)/(K4))^2 + (C(8)/(K5))^2)) - k2*C(3);

% P97 protein
dC(4) = p2*C(3) - d2D*(1-ThetaPhyA)*C(4) - d2L*ThetaPhyA*C(4);

% P51 mRNA
dC(5) =v3/(1 + (C(2)/(K6))^2 + (C(6)/(K7))^2) - k3*C(5);

% P51 protein
dC(6) =( p3*C(5) - d3D*(1-ThetaPhyA)*C(6) - d3L*ThetaPhyA*C(6));

% EL mRNA
dC(7) =(v4*ThetaPhyA/(1 + (C(2)/(K8))^2 + (C(6)/(K9))^2 + (C(8)/(K10))^2) - k4*C(7));
% dC(7) = v4*CombL/(1 + (C(2)/(K8))^2 + (C(6)/(K9))^2 + (C(8)/(K10))^2) - k4*C(7);

% EL protein
dC(8) = (p4*C(7) - (de1+(de2*C(15)+de3*C(16)+de4*C(17)+de5*C(18))/(C(15)+C(16)+C(17)+C(18)))*C(8));

% PHY A
dC(9) = (1 - ThetaPhyA)*Ap3 - (Am7*C(9)/(Ak7 + C(9))) - q2*ThetaPhyA*C(9)-kmpac*ThetaPhyA*C(9)*C(15)+kd*(1-ThetaPhyA)*C(16);

% PIF mRNA
dC(10) = v5/(1 + (C(8)/(K11))^2) - k5*C(10);

% PIF protein
dC(11) = p5*C(10) - d5D*(1-ThetaPhyA)*C(11) - d5L*ThetaPhyA*C(11);

% HYP protein
dC(12) = g1 + (g2*C(11)^2)/(K12^2 + C(11)^2);

% PHY B
dC(13) = Bp4 - ((Bm8*C(13))/(Bk8 + C(13)))-kmpbc*ThetaPhyB*C(13)*C(15)+kd*(1-ThetaPhyB)*C(17);

% CRY1
dC(14) = Cp5 - ((Cm9*C(14))/(Ck9 + C(14)))-kmcc*ThetaCry1*C(14)*C(15)+kd*(1-ThetaCry1)*C(18);

%COP1
dC(15) = -kmpac*ThetaPhyA*C(9)*C(15)+kd*(1-ThetaPhyA)*C(16)-kmpbc*ThetaPhyB*C(13)*C(15)+kd*(1-ThetaPhyB)*C(17)-kmcc*ThetaCry1*C(14)*C(15)+kd*(1-ThetaCry1)*C(18)+ (Am7*C(16)/(Ak7 + C(16))) + q2*ThetaPhyA*C(16)+((Bm8*C(17))/(Bk8 + C(17)))+((Cm9*C(18))/(Ck9 + C(18)));

%COP1:PHYA
dC(16) = kmpac*ThetaPhyA*C(9)*C(15)-kd*(1-ThetaPhyA)*C(16)- (Am7*C(16)/(Ak7 + C(16)))- q2*ThetaPhyA*C(16);

%COP1:PHYB
dC(17) = kmpbc*ThetaPhyB*C(13)*C(15)-kd*(1-ThetaPhyB)*C(17)- ((Bm8*C(17))/(Bk8 + C(17)));

%COP1:CRY1
dC(18) = kmcc*ThetaCry1*C(14)*C(15)-kd*(1-ThetaCry1)*C(18)- ((Cm9*C(18))/(Ck9 + C(18)));