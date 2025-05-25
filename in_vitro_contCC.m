close all
clear all

%% load data
% set a folder with this script inside and all the cont. CC recordings as
% .mat files that you want to analyse

files = dir ('in_vitro_contCC_example.mat');

% ToPlot = 'yes';
GroupName = 'contCC';

nCells =size(files);
nCells (2) = [];

for i=1:nCells
    temp = load(files(i).name);
    name = files(i).name;
    contCC(i).cell.name = name(1:end-4);
    contCC(i).cell.trace = temp.(char(fieldnames (temp)));
    contCC(i).cell.time = contCC(i).cell.trace (:,1);
    contCC(i).cell.voltage = contCC(i).cell.trace (:,2);
    % contCC(i).cell.stim = contCC(i).cell.trace (:,3);
end

clear files i temp name


% Set folders
mkdir(['Stat_' char(GroupName)]);                   % creates a new directory with for the analysis
oldFolder = cd(['Stat_' char(GroupName)]);          % Sets the return to old directory
addpath(pwd);                                       % adds the current path to matlab's search folders


%% analysis

results= zeros(nCells,8);

for k = 1:nCells
	% Create a mat filename, and load it into a structure called matData.
    % recording = load(sprintf('191016p%d.mat', k));
    % name_of_the_trace = fieldnames(recording);
    trace = contCC(k).cell.trace;

trace_length = 10;

Fs = length(trace)./trace_length;
ci = 5; %confident interval in ms

%
%Transforming the trace into the first derivative
time= trace(:,1);
voltage = trace(:,2);
first_derivative = diff(voltage)./diff(time); 

contCC(k).cell.deriv1 = first_derivative;

if length(first_derivative) ~= length(time)
    while length(first_derivative) ~= length(time)
        first_derivative(length(first_derivative)+1,1) = 0;
    end
end

Result.(sprintf('Cell_%d',k)).First_derivative = first_derivative;
Result.(sprintf('Cell_%d',k)).Trace = trace;


f = figure;
subplot (2,2,1)
hold on
xlim ([0 10]);
% ylim ([-50 50]);
plot(trace(:,1),first_derivative(:,1));
% fname = [contCC(k).cell.name '_deriv_vs_time'];
% saveas (f, fname, 'jpg');
hold off

subplot (2,2,2)
hold on
plot(trace(:,2),first_derivative(:,1));
% fname = [contCC(k).cell.name '_facePlot'];
% saveas (f,fname, 'jpg');
hold off

%%
% Finding thresholds in the first derivative
timestamp = 1;
threshold_times= [];
number_of_thresholds = 0;

while timestamp < (Fs*trace_length-round(Fs*ci./1000))
    mVms = first_derivative(timestamp,1);
    if mVms >= 10 && max(trace(timestamp:timestamp+round(Fs*ci./1000),2))>0
        number_of_thresholds = number_of_thresholds + 1;
        threshold_times(number_of_thresholds) = timestamp;
        timestamp = timestamp + round(Fs*ci/1000);
    else
        timestamp = timestamp + 1;
    end
end

%Firing_frequency = number_of_thresholds./trace_length

%%
%Finding the maxima and minima in the first derivative
max_first_derivative = [];
min_first_derivative = [];

for m = 1:length(threshold_times)
   tb = threshold_times(m);
   te = threshold_times(m) + round(ci*Fs/1000);
   max_first_derivative(m) = max(first_derivative(tb:te));
   min_first_derivative(m) = min(first_derivative(tb:te));
end

% max_dVdt = mean(max_first_derivative)
% min_dVdt = mean(min_first_derivative)

%%
%Finding voltages of the correspononding thresholds
voltages_at_threshold = [];

for i = 1:length(threshold_times)
    voltages_at_threshold(i) = trace(threshold_times(i),2); 
end

% Threshold  = mean(voltages_at_threshold)*1000

%%
%Measuring spike width at the threshold
all_spike_widths = [];

for j = 1:length(threshold_times)
    threshold_of_the_current_spike = voltages_at_threshold(j);
    time_of_the_current_spike = threshold_times(j);
    corresponding_voltage = 0;
    t = time_of_the_current_spike+1;
    
    while corresponding_voltage > threshold_of_the_current_spike && t < length(trace)
        corresponding_voltage = trace(t,2);
        if corresponding_voltage < threshold_of_the_current_spike
            time_of_the_corresponding_voltage = t;
            spike_width_of_the_current_spike = time_of_the_corresponding_voltage - time_of_the_current_spike;
            all_spike_widths(j) = spike_width_of_the_current_spike./Fs;
        else t = t+1;
        end  
    end
end

%Spike_width = mean(all_spike_widths)*1000 %Mean spike width in ms
  
%%
%Calculating CV
ISI = [];

for isi = 1: length(threshold_times)-1
    ISI(isi) = (threshold_times(isi+1)-(threshold_times(isi)))./Fs;
end

%CV = std(ISI)./mean(ISI)*100


%%
%Finding peaks
peak_times_in_correlation_to_threshold= [];
peak_times= [];

for i = 1:length(threshold_times)-1  %The last peak will be excluded because the AHP might not be in the trace
    peak = max(trace(threshold_times(i):threshold_times(i+1),2));
    peak_times_in_correlation_to_threshold(i) = mean(find(trace(threshold_times(i):threshold_times(i+1),2)==peak));
    peak_times(i) = round(threshold_times(i)+peak_times_in_correlation_to_threshold(i));      
end

%%
%Finding the times of the minAHP
minAHP_times_in_correlation_to_threshold= [];
minAHP_times= [];
AHP_values = [];

for i = 1:length(threshold_times) -1 %The last AHP will be excluded because it might not be in the trace
    AHP = min(trace(threshold_times(i):threshold_times(i+1),2));
    AHP_values(i) = AHP;
    minAHP_times_in_correlation_to_threshold(i) = mean(find(trace(threshold_times(i):threshold_times(i+1),2)==AHP));
    minAHP_times(i) = round(threshold_times(i)+minAHP_times_in_correlation_to_threshold(i));      
end

%minAHP = mean(AHP_values)*1000


%%
%Time from peak to minAHP
times_from_peak_to_minAHP = [];
for t = 1:length(peak_times)
    times_from_peak_to_minAHP(t) = (minAHP_times(t) - peak_times(t));
end
times_from_peak_to_minAHP;
Peak_to_minAHP= mean(times_from_peak_to_minAHP)/Fs*1000; %in ms

Firing_frequency = number_of_thresholds./trace_length;
Threshold  = mean(voltages_at_threshold)*1000;
Spike_width = mean(all_spike_widths)*1000;
meanAHP = mean(AHP_values)*1000;
max_dVdt = mean(max_first_derivative);
min_dVdt = mean(min_first_derivative);
CV = std(ISI)./mean(ISI)*100;
meanVolt = mean(voltage);


% f = figure
subplot (2,2,[3,4])
hold on
xlim ([0 10]);
ylim([-0.1 0.1]);
plot (trace(:,1), trace(:,2));
line ([0 10], [meanVolt meanVolt], 'Color',[0.4 0.7 0.2]);
scatter (threshold_times/Fs, voltages_at_threshold);
scatter (minAHP_times/Fs, AHP_values);
scatter (peak_times/Fs, trace(peak_times,2));

% fname = [contCC(k).cell.name '_trace_threshold_AHP'];
title ( contCC(k).cell.name, 'Interpreter','none');
fname = [contCC(k).cell.name '_contCC_primAn'];
saveas (f, fname, 'fig')
saveas (f, fname, 'emf')
saveas (f, fname, 'svg')
hold off
close (f)


results(k,1)= Firing_frequency;
results(k,2)= Threshold;
results(k,3)= Spike_width;
results(k,4)= meanAHP;
results(k,5)= max_dVdt;
results(k,6)= min_dVdt;
results(k,7)= CV;
results(k,8)= meanVolt;

contCC(k).cell.peak_times = peak_times/Fs;
contCC(k).cell.trace_peak_times = trace(peak_times,2);
contCC(k).cell.ISIs = ISI;
contCC(k).cell.Freq_mean = mean(1./ISI);
contCC(k).cell.CV = CV;
contCC(k).cell.Threshold_MEAN = Threshold;
contCC(k).cell.Thresholds = voltages_at_threshold;
contCC(k).cell.Spike_width_MEAN = Spike_width;
contCC(k).cell.Spike_widths = all_spike_widths;
contCC(k).cell.meanAHP = meanAHP;
contCC(k).cell.AHPs = AHP_values;
contCC(k).cell.max_dVdt = max_dVdt;
contCC(k).cell.min_dVdt = min_dVdt;
contCC(k).cell.meanVolt = meanVolt*1000;
contCC(k).cell.Spike_amplitude_mean = mean(trace(peak_times,2))*1000;
contCC(k).cell.first_derivative = first_derivative;


end

% results;

%% RASTER Plot

for j=1:nCells
    TS = contCC(j).cell.peak_times;
    name = contCC(j).cell.name;
    % Set parameters
    timeslice = 10;
    timeStartPoint = floor(min(TS));
    timeEndPoint = floor(max(TS)) + 1;  
    timetotal = timeEndPoint - timeStartPoint;
    % ls = (timetotal/timeslice);
    ls = 1;
    
    f = figure;
        set(f, 'Position', [50 50 800 100]);
        set(gca,'FontName','calibri','FontSize',12, 'YDir','reverse', 'TickDir', 'out');
        hold on
            t1 = 0;
            t2 = timeslice;
            
            for i=1:length(TS)
                line([t2+TS(i)-timetotal-timeStartPoint t2+TS(i)-timetotal-timeStartPoint],[ls-1-0.42+timeStartPoint/timeslice ls-1+0.42+timeStartPoint/timeslice],'Color',[0 0 0], 'LineWidth',0.017);% RPlot =  %
            end
        
       
            ylim([timeStartPoint/timeslice-0.9 timeStartPoint/timeslice+ls-0.1]);
            xlim([-0.1 10.1]);
            xlabel('Time [s]'); % Time is in millisecond
            set(gca,'YTickLabel',[]);
        hold off
        fname = [name '_PM_raster'];
        saveas (f, fname, 'fig')
        saveas (f, fname, 'svg')
    close (f)
end

%% ISI histogram for this group

All_ISIs = [];
for j=1:nCells
    All_ISIs = [All_ISIs contCC(j).cell.ISIs];
end

f = figure;
set (f, 'Position', [50 100 800 600]);
hold on 
xlim ([0 1.5]); xlabel ('ISI (s)');
ylim([0 150]); ylabel ('#');
histogram  (All_ISIs, 50);
hold off
fname = [GroupName '_ISI_histo'];
saveas (f, fname, 'fig')
saveas (f, fname, 'svg')
close (f)


%% Trace vs Raster 

f = figure;
set (f, 'Position', [50 50 1000 nCells*100]);

n = 1;
for j = 1:nCells
    subplot (nCells,2, n)
    hold on
        xlim ([0 10]);
        ylim([-0.1 0.05]);
        line ([0 10], [0 0], 'Color',[0.7 0.7 0.7]);
        line ([0 10], [contCC(j).cell.meanVolt contCC(j).cell.meanVolt], 'Color',[0.4 0.7 0.2]);
        plot (contCC(j).cell.trace(:,1), contCC(j).cell.trace(:,2));
%         set(gca,'XTickLabel',[]);
%         set(gca,'YTickLabel',[]);
        set(gca,'FontSize', 5);
        title (contCC(j).cell.name, 'Interpreter','none', 'FontSize', 7);
    hold off
    n = n+1;
    subplot (nCells,2, n)
    hold on
        TS = contCC(j).cell.peak_times;
        timeslice = 10;
        timeStartPoint = floor(min(TS));
        timeEndPoint = floor(max(TS)) + 1;  
        timetotal = timeEndPoint - timeStartPoint;
        ls = 1;
        t1 = 0;
        t2 = timeslice;
        for i=1:length(TS)
            line([t2+TS(i)-timetotal-timeStartPoint t2+TS(i)-timetotal-timeStartPoint],[ls-1-0.42+timeStartPoint/timeslice ls-1+0.42+timeStartPoint/timeslice],'Color',[0 0 0], 'LineWidth',0.017);% RPlot =  %
        end
        ylim([timeStartPoint/timeslice-0.9 timeStartPoint/timeslice+ls-0.1]);
        xlim([-0.1 10.1]);
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);
        set(gca,'FontSize', 5);
        % xlabel('Time [s]'); % Time is in millisecond
    hold off
    n=n+1;
end

fname = [GroupName '_ALL_PM_raster'];
saveas (f, fname, 'fig')
saveas (f, fname, 'svg')


%% trace vs Face plot



f = figure;
set (f, 'Position', [10 10 1000 nCells*100] );

n = 1;
for j = 1:nCells
    subplot (nCells,2, n)
    hold on
        xlim ([0 10]);
        ylim([-0.1 0.05]);
        line ([0 10], [0 0], 'Color',[0.7 0.7 0.7]);
        line ([0 10], [contCC(j).cell.meanVolt contCC(j).cell.meanVolt], 'Color',[0.4 0.7 0.2]);
        plot (contCC(j).cell.trace(:,1), contCC(j).cell.trace(:,2));
        set(gca,'FontSize', 5);
        title (contCC(j).cell.name, 'Interpreter','none', 'FontSize', 5);
    hold off
    n = n+1;
    subplot (nCells,2, n)
    hold on
        plot (contCC(j).cell.trace(:,2), contCC(j).cell.first_derivative(:,1));
        ylim ([-50 100]);
        xlim([-0.1 0.05]);
        set(gca,'FontSize', 5);
    hold off
    n=n+1;
end

fname = [GroupName '_ALL_PM_FACEplot'];
saveas (f, fname, 'fig')
saveas (f, fname, 'svg')


%% export data in excel file

tabl = zeros(nCells,10);

for i = 1:nCells
    tabl (i+1, 2) = contCC(i).cell.Freq_mean;
    tabl (i+1, 3) = contCC(i).cell.meanVolt;
    tabl (i+1, 4) = contCC(i).cell.CV;
    tabl (i+1, 5) = contCC(i).cell.meanAHP;
    tabl (i+1, 6) = contCC(i).cell.Threshold_MEAN;
    tabl (i+1, 7) = contCC(i).cell.max_dVdt;
    tabl (i+1, 8) = contCC(i).cell.min_dVdt;
    tabl (i+1, 9) = contCC(i).cell.Spike_width_MEAN;
    tabl (i+1, 10) = contCC(i).cell.Spike_amplitude_mean;
end

filename = [GroupName '_contCC_groupAnalysis.xlsx'];
xlswrite (filename, tabl);

tabl4CellNames = repmat (blanks(12),nCells,1);
for i = 1:nCells
    NameLength = length (contCC(i).cell.name);
    tabl4CellNames (i, 1:NameLength) = contCC(i).cell.name;
end

fileID = fopen('CellName.txt','wt');
for r=1:size(tabl4CellNames,1)
    fprintf(fileID,'%s\n',tabl4CellNames(r,:));
end
fclose(fileID);

save ([GroupName '_contCC_groupAnalysis.mat'])
close all;
