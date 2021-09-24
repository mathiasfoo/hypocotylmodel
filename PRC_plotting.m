clc
clear all
close all

global Ired Iblue

timestep=0.1;

C = 1*ones(1,18); %Initial condition for clock genes
C(12)=0;          %Initial condition for HYP - Hypocotyl Length
C(16)=0;          %Initial condition for COP1:PhyA
C(17)=0;          %Initial condition for COP1:PhyB
C(18)=0;          %Initial condition for COP1:Cry1
t1=0;             %Simulation Start Time
t2=775;           %Simulation End Time
tspan = t1:timestep:t2;

%%Background Light
Ired = 0;
Iblue =0;
[T,C] = ode15s('proposed_RBLight_ODEmod_com_full',tspan,C(end,:));

%% Free running period 
unperturbedLHY=C(:,1); % Unpertrubed LHY mRNA trajectory
[pks,locs] = findpeaks(unperturbedLHY,'MinPeakProminence',0.01); %find maxima points 
peaktime=T(locs);
freeDD=peaktime(end)-peaktime(end-1);

%% Set the reference peak time of the circadian gene expression before light stimulus
findstartpoint=peaktime>200; %refence point right before light stimuls
startpoint=min(peaktime(findstartpoint));
findreference=peaktime>682; %reference point for non-stimulated circadian rhythm 
referpoint=min(peaktime(findreference));

%%Light stimulus
stduration =1; % stimulus duration
phaseshift=[];
phasechange=[];
for i =startpoint+[0:2:freeDD]
     
    %Initial conditions
    C0 = 1*ones(1,18);
    C0(12)=0;
    C0(16)=0;
    C0(17)=0;
    C0(18)=0;
    
    % before stimulus
    t1=0;
    t2=i;
    tspan = t1:timestep:t2;
    Ired = 0;
    Iblue = 0;   
    [T1,C1] = ode15s('proposed_RBLight_ODEmod_com_full',tspan,C0(end,:));
    
    % during stimulus
    t1=i;
    t2=i+stduration;
    tspan = t1:timestep:t2;
    Ired = 40;
    Iblue= 0;   
    [T2,C2] = ode15s('proposed_RBLight_ODEmod_com_full',tspan,C1(end,:));
    
    % after stimulus
    t1=i+stduration;
    t2=775;
    tspan = t1:timestep:t2;    
    Ired = 0;
    Iblue = 0;   
    [T3,C3] = ode15s('proposed_RBLight_ODEmod_com_full',tspan,C2(end,:));

    T=[T1;T2;T3];
    C=[C1;C2;C3];

    [pks,locs]=findpeaks(C(:,1),'MinPeakProminence',0.01);
    ppeaktimes=T(locs);
    
    % Estimation of phaseshift and phase of stimulus
    findppeaktime=abs(referpoint - ppeaktimes) < freeDD/2;
    ppeaktime=(ppeaktimes(findppeaktime));
    phaseshift=[phaseshift,referpoint-ppeaktime];
    phasechange=[phasechange,i-startpoint]; %phase of stimulus
 
end

%%Plot PRC
figure('Position', [800 400 250 200])
plot(phasechange/freeDD,phaseshift/freeDD,'k','LineWidth',6)
hold on
plot(0:0.1:1,zeros(1,11),'k-','LineWidth',2)
xticks([0 0.25 0.5 0.75 1])
yticks([-0.6,-0.4,-0.2,0,0.2,0.4,0.6]);
axis([0 1 -0.6 0.6]);
ylabel('Phase Shift (rad/2\pi)');
xlabel('Normalised Phase of Stimulus');

