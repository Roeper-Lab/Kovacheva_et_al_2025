clear all
close all

load ('behavior_GroupsRotAnalysis.mat');


allRelRotRi = absNumRi*100./rotAbsNum;
OHDAall = allRelRotRi(:,13:22);
ACSFall = allRelRotRi(:,1:12);


sessions = [-12 -8 -4 4 8 12 16 20 24 28 32 36 40 44 48 52 56 60 64 68]; %  days of experiment
p50 = ones (20); % p50 is for drawing a line at 50%
p50 = p50*50;

figure;
hold on;

xlim ([-12 64]); %limiting axes length
ylim ([0 100]);
line([-12 64], [50 50], 'Color', 'k', 'LineStyle', ':'); 
xlabel ('Sessions (Days)');
ylabel ('%');
title('Percent of right rotations');

meanR  = mean (OHDAall'); %the mean for 6-OHDA mice % right rotations
std_OHDAall = std(OHDAall'); % the standard error for 6-OHDA mice % right rotations
SEM_OHDAall = std_OHDAall / sqrt (length(OHDAall)); % calculating the SEM for 6-OHDA mice
SEM_OHDAall = SEM_OHDAall*1.96;


plot(sessions, meanR, 'Color', 'r')

endV = length(OHDAall) - 1; 
for i = 1:endV % shading the SEM for 6-OHDA
    x2(1) = sessions (i);
    y2(1) = meanR (i) - SEM_OHDAall(i);
    x2(2) = sessions (i);
    y2(2) = meanR (i) + SEM_OHDAall(i);
    x2(3) = sessions (i+1);
    y2(3) = meanR (i+1) + SEM_OHDAall(i+1);
    x2(4) = sessions (i+1);
    y2(4) = meanR (i+1) - SEM_OHDAall(i+1);
    patch (x2, y2,'red', 'FaceAlpha', .2, 'LineStyle', 'none')
end
 
meanA  = mean (ACSFall'); % all the same as 6-OHDA mice but for the ACSF mice
std_ACSFall = std(ACSFall');
SEM_ACSFall = std_ACSFall / sqrt (length(ACSFall));
SEM_ACSFall = SEM_ACSFall*1.96;
sessions = [-12 -8 -4 4 8 12 16 20 24 28 32 36 40 44 48 52 56 60 64 68];

plot(sessions, meanA, 'blue')

endA = length(ACSFall) - 1;
for i = 1:endA
    x(1) = sessions (i);
    y(1) = meanA (i) - SEM_ACSFall(i);
    x(2) = sessions (i);
    y(2) = meanA (i) + SEM_ACSFall(i);
    x(3) = sessions (i+1);
    y(3) = meanA (i+1) + SEM_ACSFall(i+1);
    x(4) = sessions (i+1);
    y(4) = meanA (i+1) - SEM_ACSFall(i+1);
    patch (x, y,'blue', 'FaceAlpha', .2, 'LineStyle', 'none')
end
hold off;

% save ('PercentRightRotations_Analysis_plot.mat');