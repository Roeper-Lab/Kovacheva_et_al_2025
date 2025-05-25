clear all
close all

%% load data

files = dir ('in_vito_onCell_example.mat');

% ToPlot = 'yes';
GroupName = 'onCell';

nCells =size(files);
nCells (2) = [];

for i=1:nCells
    temp = load(files(i).name);
    name = files(i).name;
    onCell(i).cell.name = name(1:end-4);
    onCell(i).cell.trace = temp.(char(fieldnames (temp)));
    onCell(i).cell.time = onCell(i).cell.trace (:,1);
    onCell(i).cell.voltage = onCell(i).cell.trace (:,2);
end

clear files i temp name


% Set folders
mkdir(['Stat_' char(GroupName)]);                   % creates a new directory with for the analysis
oldFolder = cd(['Stat_' char(GroupName)]);          % Sets the return to old directory
addpath(pwd);                                       % adds the current path to matlab's search folders

%% analysis in an INVERTED trace!!!!! ---> all mV-values  inverted!!!!

for i = 1:nCells
    % low pass filter parameters 
    ci = 5; %confident interval in ms
    fs = length(onCell(i).cell.trace(:,1))./onCell(i).cell.trace(end,1); % sampling rate 
    HFcut = 50; % high-cut frequency for lowpass filter
    order = 2; % for highpass filter
    
    % high PASS FILTER - to remove high frequency noise
    [b,a] = butter(order,HFcut/(fs/2),'high');
    CorrTrace = filtfilt(b,a, onCell(i).cell.trace(:,2));
    
    % CorrTrace = onCell(i).cell.trace(:,2) - mean (onCell(i).cell.trace(:,2)); % substract offset (better visualization in plot)
    [pks,locs] = findpeaks(-1.*CorrTrace);  % trace inverted!!!!!
    maxPeak = max(-1.*CorrTrace);
    thresh = 0.4*maxPeak;
    spikes_mV = pks(pks>thresh);
    spikes_s = onCell(i).cell.trace(locs(pks>thresh),1);

    for j = 1:length(spikes_mV)-1
        if j > length(spikes_mV)-1 
            break;
        end
        if 1/(spikes_s(j+1)-spikes_s(j)) > 100 % find two points at one spike
            if spikes_mV(j+1) > spikes_mV(j)
                spikes_s(j)= [];
                spikes_mV(j)= [];
                if j > length(spikes_mV)-1 
                    break;
                end
            else
                spikes_s(j+1)= [];
                spikes_mV(j+1)= [];
                if j > length(spikes_mV)-1 
                    break;
                end
            end
        end
    end
    
        
    f = figure ('Position', [50, 50, 1200, 600]);
    subplot(2,1,1)
    hold on
        plot(onCell(i).cell.trace(:,1), onCell(i).cell.trace(:,2));
        title([onCell(i).cell.name ' - original'], 'interpreter','none')
    hold off
    subplot(2,1,2)
    hold on
        scatter (spikes_s, -1.*spikes_mV);
        plot(onCell(i).cell.trace(:,1), CorrTrace);
        title([onCell(i).cell.name ' - high pass filter'],'interpreter','none')
    hold off
    fname = [onCell(i).cell.name '_onCell'];
    saveas (f, fname, 'fig')
    saveas (f, fname, 'svg')
    saveas (f, fname, 'jpg')
    
    
    ISIs = diff(spikes_s);
    freq = 1./ISIs;
    meanFreq = mean(freq);
    CV = std(ISIs)./mean(ISIs)*100;
    
    onCell(i).cell.thresh = thresh;
    onCell(i).cell.spikes_mV = spikes_mV;
    onCell(i).cell.spikes_s = spikes_s;
    onCell(i).cell.ISIs = ISIs;
    onCell(i).cell.freq = freq;
    onCell(i).cell.meanFreq = meanFreq;
    onCell(i).cell.CV = CV;
    
end

% close all

%% export data in excel file

tabl = zeros(nCells,10);

for i = 1:nCells
    tabl (i+1, 2) = onCell(i).cell.meanFreq;
    tabl (i+1, 3) = onCell(i).cell.CV;
end

filename = [GroupName '_onCell_groupAnalysis.xlsx'];
xlswrite (filename, tabl);

tabl4CellNames = repmat (blanks(12),nCells,1);
for i = 1:nCells
    NameLength = length (onCell(i).cell.name);
    tabl4CellNames (i, 1:NameLength) = onCell(i).cell.name;
end

fileID = fopen('CellName.txt','wt');
for r=1:size(tabl4CellNames,1)
    fprintf(fileID,'%s\n',tabl4CellNames(r,:));
end
fclose(fileID);

save ([GroupName '_onCell_groupAnalysis.mat']);
close all;