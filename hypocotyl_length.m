clc
clear all
close all

global Ired Iblue eta1 eta2 Ctot
pPlot = [];
Nday = 10; % Experimental data collected at Day 10

for LH=[18:2:22] %Photoperiod (hours in darkness)
    for wave = 1:1:3 %light colours
    % Set the Intensity to 26.62 for ON and 0 for OFF
    
    if wave == 1
    %Blue light
    IntensityBB = 1*26.62*[ones(1,24-LH) zeros(1,LH)];
    IntensityRR = 1*0*[ones(1,24-LH) zeros(1,LH)];
    
    elseif wave == 2
    %Red light 
    IntensityBB = 1*0*[ones(1,24-LH) zeros(1,LH)];
    IntensityRR = 1*26.62*[ones(1,24-LH) zeros(1,LH)]; 
    
    else 
    %Mixed
    IntensityBB = 1*26.62*[ones(1,24-LH) zeros(1,LH)];
    IntensityRR = 1*26.62*[ones(1,24-LH) zeros(1,LH)];
    end
    
    IntensityBB = repmat(IntensityBB,1,Nday);
    IntensityRR = repmat(IntensityRR,1,Nday);
    
    ProteinLevel = [];
    C = 1*ones(1,18);
    C(12)=0;
    C(16)=0;
    C(17)=0;
    C(18)=0;
    for t = 1:length(IntensityRR)
        tspan = [t t+1];
        Ired = IntensityRR(t);
        Iblue = IntensityBB(t);
        [T,C] = ode15s('proposed_RBLight_ODEmod_com_full',tspan,C(end,:));
        ProteinLevel = [ProteinLevel; C(end,:)];
    end
    
    % Calculating hypocotyl length
    tp = 1:length(IntensityRR);
    hyplength(LH) = ProteinLevel(10*24,12); % Hypocotyl length on Day 10
    hyplength(hyplength== 0) = [];  
    end
end

%% Simulated Hypocotyl Length Plotting
h=round(hyplength,2); %round off simulated hypocotyl length to 2 decimal places
S_hyp=[h(1:3);h(4:6);h(7:9)];
figure
b = bar(S_hyp, 'grouped');
ylim([0 20]);
set(gca,'XTickLabel',{'6L18D','4L20D','2L22D'})
legend('Blue','Red','Blue+Red')
ylabel('Hypocotyl Length (mm)');
xlabel('Light Duration (h)');
title('Simulated Hypocotyl Length');

%%Measured Hypocotyl Length Plotting
M_hyp=[5.13 7.01 5.59;6.14 9.42 6.95;8.77 10.6 9.2];
sd_hyp=[0.87 2.05 1;1.14 1.72 1.19;2.60 2.18 1.57];
figure
b = bar(M_hyp, 'grouped');
hold on
ngroups = size(M_hyp, 1);
nbars = size(M_hyp, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, M_hyp(:,i), sd_hyp(:,i), '.');
end
hold off
ylim([0 20]);
ylabel('Hypocotyl Length (mm)');
xlabel('Light Duration (h)');
set(gca,'XTickLabel',{'6L18D','4L20D','2L22D'})
legend('Blue','Red','Blue+Red')
title('Measured Hypocotyl Length');

%%Min-max Scaling Simulated Hypocotyl Length Plotting 
LH1=([8.21 8.86  8.50]-8.21).*2.8923+5.13;
LH2=([11.98 12.67 12.27]-11.98).*4.7536+6.14;
LH3=([ 15.88 16.14 16.01]-15.88).*7.0385+8.77;
N_hyp=[LH1;LH2;LH3];
figure
b = bar(N_hyp, 'grouped');
ylim([0 20]);
set(gca,'XTickLabel',{'6L18D','4L20D','2L22D'})
legend('Blue','Red','Blue+Red')
ylabel('Hypocotyl Length (mm)');
xlabel('Light Duration (h)');
title('Simulated Hypocotyl Length (Normalised)');





