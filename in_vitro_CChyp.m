close all
clear all

files = dir ('in_vitro_CChyp_example.mat');

ToPlot = 'yes';
GroupName = 'CChyp';

nCells = size(files);
nCells (2) = [];

for i=1:nCells
    temp = load(files(i).name);
    name = files(i).name;
    HypCC(i).cell.name = name(1:end-4);
    tracesNames = fieldnames (temp); % acquiring the traces names
    HypCC(i).cell.tracesNames = tracesNames;
    nTraces = size (HypCC(i).cell.tracesNames);
    nTraces (2) = [];
    HypCC(i).cell.nTraces = nTraces;
    for j=1:nTraces
        HypCC(i).cell.(char(tracesNames(j,:)))= temp.(char(tracesNames(j, :)));
    end
end

clearvars -except GroupName nCells HypCC ToPlot 

% Set folders
mkdir(['Stat_' char(GroupName)]);                   % creates a new directory with for the analysis
oldFolder = cd(['Stat_' char(GroupName)]);          % Sets the return to old directory
addpath(pwd); 

%% 
number_of_cells = nCells;


for k = 1:nCells %1:number_of_cells %LOOP FOR IMPORTING CELLS
    % recording = load(sprintf('cell_%d.mat', k));
    % name_of_the_trace = fieldnames(recording);
    recording = HypCC(k).cell; % load(sprintf('191016p%d.mat', k));
    name_of_the_trace = HypCC(k).cell.tracesNames;
    T = 0;
    
    T_n = 1;
    
    for r = 1:length(name_of_the_trace) %LOOP FOR IMPORTING SINGLE TRACES
       T = T + 1;
       trace = recording(1).(char(name_of_the_trace(r)));
       stimulus(:,1) = trace(:,1);
       stimulus(:,2) = trace(:,3);
       
       if strcmp (ToPlot, 'yes')
            f = figure; %Plotting trace & stimulus
            hold on;
            subplot(3,1,1);
            plot(trace(:,1),trace(:,2));
            title(['Recording trace ' char(recording.name) ' ' char(name_of_the_trace(r))], 'Interpreter', 'none');
            subplot(3,1,2);
            plot(stimulus(:,1),stimulus(:,2));
            title('Stimulus trace');
       end

    Fs = 20000;
    tl = length(trace)/Fs; %trace length in s
    ci = 5*Fs/1000;        %confident interval 5ms?

    if max(abs(stimulus(:,2))) > 1e-12 %if stimulus is bigger than 1pA -> trace will be analyzed

    sb = find(abs(stimulus(:,2)) > 1e-12 ==1,1,'first');
    se = find(abs(stimulus(:,2)) > 1e-12 ==1,1,'last');
    
        if strcmp(ToPlot, 'yes')
            subplot(3,1,2);
            hold on;
            plot(sb/Fs,0,'ro');
            plot(se/Fs,0,'ro');
        end

    fd = [diff(trace(:,2))./diff(trace(:,1));0];
    trace(:,3) = 1:length(trace); %filling the timestamps into the 3rd row 
    
    if strcmp (ToPlot, 'yes')
        subplot (3,1,3)
        hold on;
        plot(trace(:,1),fd);
        title('First derivative trace');
    end
    
    %% Finding thresholds

    fdo10 = trace(fd>=20,1:3); %finding all values where fd is >=20 !!!!!
    
    if isempty(fdo10)
        continue;
    end
    
    fdo10 = fdo10([1;find(diff(fdo10(:,1))>=0.003)+1],:);

    tn = 0; %number of thresholds
    for i=1:length(fdo10(:,1))
          if  fdo10(i,3)+ci*4 < length(trace) && ...
              max(trace((fdo10(i,3)):(fdo10(i,3)+ci),2))> -0.02 && ...
              min(fd((fdo10(i,3)):(fdo10(i,3)+ci),1))< -5  
          tn = tn + 1;
          thresholds(tn,1:3) = fdo10(i,1:3);
          end    
    end


    
    %% Calculating other values 
    [ttemp, ~] = size(thresholds);
    for i = 1:ttemp % length(thresholds)-1
        if thresholds(i,3)+ci*40 > 200000
            break
        else
            [~,x] = max(trace((thresholds(i,3)):(thresholds(i,3)+ci),2));
            peaks(i,1:3) = trace(thresholds(i,3)+x,1:3);
            [~,x] = min(trace((thresholds(i,3)+ci):(thresholds(i,3)+ci*40),2));
            minAHPs(i,1:3) = trace(thresholds(i,3)+ci+x,1:3);
            x = find((trace((thresholds(i,3)):(thresholds(i,3)+ci*2),2) <= thresholds(i,2))==1,1);
            SWs(i,1:3) = trace(thresholds(i,3)+x,1:3);
            spike_widths(i) = x/Fs*1000; %in ms
        end
    end

    
    ISIs = diff(peaks(:,1));


    %%  slope = (y2-y1)/(x2-x1); % x - time; y - Voltage
    
    %% SLOPE -52 till -45 mV
    
    %find point -52mV
    at52mv = [];
    for i = se:length(trace(:,1))
        % at52mV - index, time, mV at -52mV after CChyp
        if trace(i, 2) == -0.052 || trace (i,2)>-0.052
            at52mv = [i trace(i, 1) trace(i,2)];
            break
        end
    end
    
    %find point -45mV
    at45mv = [];
    at46mv = [];
    
    for i = se:length(trace(:,1))
        % at45mV - index, time, mV at -45mV after CChyp
        if trace(i, 2) == -0.045 || trace (i,2)>-0.045
            at45mv = [i trace(i,1) trace(i,2)];
            break
        end
    end
    
    if strcmp(ToPlot, 'yes')
        subplot(3,1,1);
        hold on;
        plot(at52mv(1)/Fs,at52mv(3),'Marker','x', 'Color', 'r');
        plot(at45mv(1)/Fs,at45mv(3),'Marker','x', 'Color', 'r');
        line ([at52mv(1)/Fs at45mv(1)/Fs], [at52mv(3) at45mv(3)], 'Color', 'r');
    end
    
    slope4552 = (at45mv(3) - at52mv(3))/(at45mv(2)-at52mv(2)); % mV/mS
    
    % if the slope is too long / few failed spikes ... -> check slope at -46mV
    if (at45mv(2)-at52mv(2)) > 0.3
        for i = se:length(trace(:,1))
        % at46mV - index, time, mV at -46mV after CChyp
            if trace(i, 2) == -0.046 || trace (i,2)>-0.046
                at46mv = [i trace(i,1) trace(i,2)];
                break
            end
        end
        
        if strcmp(ToPlot, 'yes')
            subplot(3,1,1);
            hold on;
            plot(at52mv(1)/Fs,at52mv(3),'Marker','x', 'Color', 'b');
            plot(at46mv(1)/Fs,at46mv(3),'Marker','x', 'Color', 'b');
            line ([at52mv(1)/Fs at46mv(1)/Fs], [at52mv(3) at46mv(3)], 'Color', 'b');
        end
        
         slope4652 = (at46mv(3) - at52mv(3))/(at46mv(2)-at52mv(2)); % mV/mS
        
    end
    
    
    
    %% Rebound delay ms
    
    PostStimSpike = [];
    for i = 1:(length(thresholds)-1)
        if thresholds(i,1)<se/Fs && thresholds(i+1,1)>se/Fs
            PostStimSpike = [thresholds(i+1,1) thresholds(i+1,2) thresholds(i+1,3)]; 
            if length(thresholds)-i>1
                Freqs_after_Hyp = 1./diff(thresholds(i+1:end, 1));
            end
            break
        elseif i==1 && thresholds(i,1)>se/Fs
            PostStimSpike = [thresholds(i,1) thresholds(i,2) thresholds(i,3)]; 
            break
        end
    end
    
    ReboundDelay = (PostStimSpike(1) - se/Fs)*1000;
    
    %% stim voltage
    
    Stim_pA = mean(stimulus(sb:se,2));
    [minSAG, IndSAG] = min(trace(sb:se,2));
    last100msStim = mean(trace((se-ci*20):se,2))*1000;
    SagV = abs( abs(last100msStim) - abs(minSAG*1000) );
    
    % tau HCN
    Sb_minSAG_time_ms = trace(sb,1) - trace(sb+IndSAG,1);
    Sb_mV = trace(sb,2);
    Sb_minSAG_mV = (abs(minSAG) - abs(Sb_mV))*1000;
    Rmemb_GOm = Sb_minSAG_mV/abs(Stim_pA); % calc C
    
    tau_HCN_mV = Sb_mV - Sb_minSAG_mV/(exp(1)*500);
    diff_temp = trace(sb:(sb+IndSAG),2) - tau_HCN_mV;
    [tau_HCN_diff_mV, tau_HCN_diff_idx] = min(abs(diff_temp));
    tau_HCN_ms = (trace(sb+tau_HCN_diff_idx,1) - trace(sb,1))*1000;
    
    
    % tau 2
    
    minSAG_se_time_ms = trace(se,1) - trace(sb+IndSAG,1);
    Se_mV = trace(se,2);
    minSAG_se_mV = (abs(minSAG) - abs(Se_mV))*1000;
    
    tau2_mV = minSAG*1000 + SagV*2/exp(1);
    diff_temp2 = trace((sb+IndSAG):se,2) - tau2_mV/1000;
    [tau2_diff_mV, tau2_diff_idx] = min(abs(diff_temp2));
    tau2_ms = (trace(sb+IndSAG+tau2_diff_idx,1) - trace(sb+IndSAG,1))*1000;
    
    if strcmp(ToPlot, 'yes')
        subplot(3,1,1);
        hold on;
        plot(thresholds(:,1),thresholds(:,2),'go');
        plot(peaks(:,1),peaks(:,2),'kx');
        plot(minAHPs(:,1),minAHPs(:,2),'m*');
        plot(SWs(:,1),SWs(:,2),'ro');
        scatter(trace(sb+IndSAG,1),minSAG); % plots min during stim
        scatter (trace(sb,1),Sb_mV); % plots begining of the stim
        scatter(trace(sb+tau_HCN_diff_idx,1), tau_HCN_mV); % plots tau-HCN, 'MarkerEdgeColor', 'k'
        scatter(trace(sb+IndSAG+tau2_diff_idx,1), tau2_mV/1000); % plots tau2, 'MarkerEdgeColor', 'y'
        fname = [HypCC(k).cell.name '_' sprintf('Trace_%d',T_n) '_CChyp'];
        saveas (f, fname, 'fig');
        saveas (f, fname, 'svg');
        close(f)
    end

    
    %% save in struct
    
    HypCC(k).cell.(sprintf('Trace_%d',T_n)).trace = trace;
    HypCC(k).cell.(sprintf('Trace_%d',T_n)).ISIs = ISIs;
    HypCC(k).cell.(sprintf('Trace_%d',T_n)).minAHPs = minAHPs;
    HypCC(k).cell.(sprintf('Trace_%d',T_n)).peaks = peaks;
    HypCC(k).cell.(sprintf('Trace_%d',T_n)).PostStimSpike = PostStimSpike;
    HypCC(k).cell.(sprintf('Trace_%d',T_n)).ReboundDelay = ReboundDelay;
    HypCC(k).cell.(sprintf('Trace_%d',T_n)).spikeWidths = spike_widths;
    HypCC(k).cell.(sprintf('Trace_%d',T_n)).SWs = SWs;
    HypCC(k).cell.(sprintf('Trace_%d',T_n)).thresholds = thresholds;
    HypCC(k).cell.(sprintf('Trace_%d',T_n)).PostStimSpikeThreshold = PostStimSpike(2)*1000; % in mV
    HypCC(k).cell.(sprintf('Trace_%d',T_n)).Stim_pA = Stim_pA*1.0e+12; % in pA
    HypCC(k).cell.(sprintf('Trace_%d',T_n)).minSAG = minSAG*1000;
    HypCC(k).cell.(sprintf('Trace_%d',T_n)).IndSAG = IndSAG;
    HypCC(k).cell.(sprintf('Trace_%d',T_n)).Sb = sb;
    HypCC(k).cell.(sprintf('Trace_%d',T_n)).Sb_minSAG_time_ms = abs(Sb_minSAG_time_ms*1000);
    HypCC(k).cell.(sprintf('Trace_%d',T_n)).Sb_mV = Sb_mV*1000;
    HypCC(k).cell.(sprintf('Trace_%d',T_n)).Sb_minSAG_mV = Sb_minSAG_mV;
    HypCC(k).cell.(sprintf('Trace_%d',T_n)).Rmemb_GOm = Rmemb_GOm;
    HypCC(k).cell.(sprintf('Trace_%d',T_n)).tau_HCN_mV = tau_HCN_mV*1000;
    HypCC(k).cell.(sprintf('Trace_%d',T_n)).tau_HCN_ms = tau_HCN_ms;
    HypCC(k).cell.(sprintf('Trace_%d',T_n)).diff_temp = diff_temp;
    HypCC(k).cell.(sprintf('Trace_%d',T_n)).tau_HCN_diff_idx = tau_HCN_diff_idx;
    HypCC(k).cell.(sprintf('Trace_%d',T_n)).tau_HCN_diff_mV = tau_HCN_diff_mV;
    HypCC(k).cell.(sprintf('Trace_%d',T_n)).last100msStim = last100msStim;
    HypCC(k).cell.(sprintf('Trace_%d',T_n)).SagV = SagV;
    HypCC(k).cell.(sprintf('Trace_%d',T_n)).minSAG_se_mV = minSAG_se_mV;
    HypCC(k).cell.(sprintf('Trace_%d',T_n)).minSAG_se_time_ms = minSAG_se_time_ms;
    HypCC(k).cell.(sprintf('Trace_%d',T_n)).diff_temp2 = diff_temp2;
    HypCC(k).cell.(sprintf('Trace_%d',T_n)).tau2_diff_idx = tau2_diff_idx;
    HypCC(k).cell.(sprintf('Trace_%d',T_n)).tau2_mV = tau2_mV;
    HypCC(k).cell.(sprintf('Trace_%d',T_n)).tau2_ms = tau2_ms;
    HypCC(k).cell.(sprintf('Trace_%d',T_n)).slope4552 = slope4552;
    HypCC(k).cell.(sprintf('Trace_%d',T_n)).at45mv = at45mv;
    HypCC(k).cell.(sprintf('Trace_%d',T_n)).at46mv = at46mv;
    HypCC(k).cell.(sprintf('Trace_%d',T_n)).at52mv = at52mv;
    HypCC(k).cell.(sprintf('Trace_%d',T_n)).sb = sb/Fs; % stimulus begin
    HypCC(k).cell.(sprintf('Trace_%d',T_n)).se = se/Fs; % stiimulus end
    if exist('slope4550', 'var')
         HypCC(k).cell.(sprintf('Trace_%d',T_n)).slope4652 = slope4652;
    end
    if exist ('Freqs_after_Hyp','var')
        HypCC(k).cell.(sprintf('Trace_%d',T_n)).Freqs_after_Hyp = Freqs_after_Hyp;
    end
    
    clearvars -except name_of_the_trace r T k recording number_of_cells HypCC T_n GroupName nCells ToPlot Fs
    
    T_n = T_n+1; 
    end
    
     
    end %LOOP FOR IMPORTING SINGLE TRACES
    HypCC(k).cell.nTraces = T_n-1;
end %LOOP FOR IMPORTING CELLS


tabl = zeros(nCells*4,16);
row = 2;
for i = 1:nCells % going through all sweeps/diferent stimulation
    
    for j = 1:4
        if isfield(HypCC(i).cell,sprintf('Trace_%d',j)) 
            tabl (row, 2) = j;
            tabl (row, 3) = HypCC(i).cell.(sprintf('Trace_%d',j)).ReboundDelay;
            tabl (row, 4) = HypCC(i).cell.(sprintf('Trace_%d',j)).Stim_pA;
            tabl (row, 5) = HypCC(i).cell.(sprintf('Trace_%d',j)).minSAG;
            tabl (row, 6) = HypCC(i).cell.(sprintf('Trace_%d',j)).last100msStim;
            tabl (row, 7) = HypCC(i).cell.(sprintf('Trace_%d',j)).SagV;
            tabl (row, 8) = HypCC(i).cell.(sprintf('Trace_%d',j)).tau_HCN_ms;
            tabl (row, 9) = HypCC(i).cell.(sprintf('Trace_%d',j)).Sb_minSAG_time_ms;
            tabl (row, 10) = HypCC(i).cell.(sprintf('Trace_%d',j)).tau_HCN_diff_mV;
            tabl (row, 11) = HypCC(i).cell.(sprintf('Trace_%d',j)).tau2_ms;
            tabl (row, 12) = HypCC(i).cell.(sprintf('Trace_%d',j)).minSAG_se_time_ms;
            tabl (row, 13) = HypCC(i).cell.(sprintf('Trace_%d',j)).slope4552;
            if isfield (HypCC(i).cell.(sprintf('Trace_%d',j)), 'slope4652')
                tabl (row, 14) = HypCC(i).cell.(sprintf('Trace_%d',j)).slope4652;
            end
            tabl (row, 15) =  HypCC(i).cell.(sprintf('Trace_%d',j)).PostStimSpikeThreshold;
            if isfield (HypCC(i).cell.(sprintf('Trace_%d',j)), 'Freqs_after_Hyp') && ~isempty (HypCC(i).cell.(sprintf('Trace_%d',j)).Freqs_after_Hyp)
                tabl (row, 16) = HypCC(i).cell.(sprintf('Trace_%d',j)).Freqs_after_Hyp(1);
            end
            row = row +1;
        end
     end
end

filename = [GroupName '_CCHyp.xlsx'];
xlswrite (filename, tabl);


% Cell names export
tabl4CellNames = repmat (blanks(12),nCells*4,1);
row = 1;
for i = 1:nCells
    for j = 1:4
        if isfield(HypCC(i).cell,sprintf('Trace_%d',j)) 
           NameLength = length (HypCC(i).cell.name);
           tabl4CellNames (row, 1:NameLength) = HypCC(i).cell.name; % 
           row = row+1;
        end        
     end
     
end

fileID = fopen('CellName.txt','wt');
for r=1:size(tabl4CellNames,1)
     fprintf(fileID,'%s\n',tabl4CellNames(r,:));
end
fclose(fileID);

% close all


%%  export for the last hypCC trace

tabl = zeros(nCells,16);
row = 2;
for i = 1:nCells % going through all sweeps/diferent stimulation
    j = HypCC(i).cell.nTraces;
    tabl (i+1, 2) = j;
    tabl (i+1, 3) = HypCC(i).cell.(sprintf('Trace_%d',j)).ReboundDelay;
    tabl (i+1, 4) = HypCC(i).cell.(sprintf('Trace_%d',j)).Stim_pA;
    tabl (i+1, 5) = HypCC(i).cell.(sprintf('Trace_%d',j)).minSAG;
    tabl (i+1, 6) = HypCC(i).cell.(sprintf('Trace_%d',j)).last100msStim;
    tabl (i+1, 7) = HypCC(i).cell.(sprintf('Trace_%d',j)).SagV;
    tabl (i+1, 8) = HypCC(i).cell.(sprintf('Trace_%d',j)).tau_HCN_ms;
    tabl (i+1, 9) = HypCC(i).cell.(sprintf('Trace_%d',j)).Sb_minSAG_time_ms;
    tabl (i+1, 10) = HypCC(i).cell.(sprintf('Trace_%d',j)).tau_HCN_diff_mV;
    tabl (i+1, 11) = HypCC(i).cell.(sprintf('Trace_%d',j)).tau2_ms;
    tabl (i+1, 12) = HypCC(i).cell.(sprintf('Trace_%d',j)).minSAG_se_time_ms;
    tabl (i+1, 13) = HypCC(i).cell.(sprintf('Trace_%d',j)).slope4552;
    
    if isfield (HypCC(i).cell.(sprintf('Trace_%d',j)), 'slope4652')
        tabl (i+1, 14) = HypCC(i).cell.(sprintf('Trace_%d',j)).slope4652;
    end
    
    tabl (i+1, 15) =  HypCC(i).cell.(sprintf('Trace_%d',j)).PostStimSpikeThreshold;
    
    if isfield (HypCC(i).cell.(sprintf('Trace_%d',j)), 'Freqs_after_Hyp') && ~isempty (HypCC(i).cell.(sprintf('Trace_%d',j)).Freqs_after_Hyp)
        tabl (i+1, 16) = HypCC(i).cell.(sprintf('Trace_%d',j)).Freqs_after_Hyp(1);
    end
end

filename = [GroupName '_CCHyp_lastTrace.xlsx'];
xlswrite (filename, tabl);

% Cell names export
tabl4CellNames = repmat (blanks(12),nCells,1);
row = 1;
for i = 1:nCells
    NameLength = length (HypCC(i).cell.name);
    tabl4CellNames (i, 1:NameLength) = HypCC(i).cell.name; % 
end
     
fileID = fopen('CellName_lastTrace.txt','wt');
for r=1:size(tabl4CellNames,1)
     fprintf(fileID,'%s\n',tabl4CellNames(r,:));
end
fclose(fileID);


%%
%% export for the closest to -80mV  & all traces CC hyp PLOT ---- 

f = figure;
set (f, 'Position', [50 50 1000 nCells*100]);

tabl = zeros(nCells,16);
row = 2;
nplot = 1;
for i = 1:nCells % going through all sweeps/diferent stimulation
    j = HypCC(i).cell.nTraces;
    minV_SAGs = [];
    if j>1
        for k=1:j
            minV_SAGs = [minV_SAGs HypCC(i).cell.(sprintf('Trace_%d',k)).minSAG];
        end
        [n, idx] = min( abs(abs(minV_SAGs) - 80) );
    else idx = 1;
    end
    
    tabl (i+1, 2) = idx;
    tabl (i+1, 3) = HypCC(i).cell.(sprintf('Trace_%d',idx)).ReboundDelay;
    tabl (i+1, 4) = HypCC(i).cell.(sprintf('Trace_%d',idx)).Stim_pA;
    tabl (i+1, 5) = HypCC(i).cell.(sprintf('Trace_%d',idx)).minSAG;
    tabl (i+1, 6) = HypCC(i).cell.(sprintf('Trace_%d',idx)).last100msStim;
    tabl (i+1, 7) = HypCC(i).cell.(sprintf('Trace_%d',idx)).SagV;
    tabl (i+1, 8) = HypCC(i).cell.(sprintf('Trace_%d',idx)).tau_HCN_ms;
    tabl (i+1, 9) = HypCC(i).cell.(sprintf('Trace_%d',idx)).Sb_minSAG_time_ms;
    tabl (i+1, 10) = HypCC(i).cell.(sprintf('Trace_%d',idx)).tau_HCN_diff_mV;
    tabl (i+1, 11) = HypCC(i).cell.(sprintf('Trace_%d',idx)).tau2_ms;
    tabl (i+1, 12) = HypCC(i).cell.(sprintf('Trace_%d',idx)).minSAG_se_time_ms;
    tabl (i+1, 13) = HypCC(i).cell.(sprintf('Trace_%d',idx)).slope4552;
    if isfield (HypCC(i).cell.(sprintf('Trace_%d',idx)), 'slope4652')
        tabl (i+1, 14) = HypCC(i).cell.(sprintf('Trace_%d',idx)).slope4652;
    end
    
    tabl (i+1, 15) =  HypCC(i).cell.(sprintf('Trace_%d',idx)).PostStimSpikeThreshold;
    
    if isfield (HypCC(i).cell.(sprintf('Trace_%d',j)), 'Freqs_after_Hyp') && ~isempty (HypCC(i).cell.(sprintf('Trace_%d',j)).Freqs_after_Hyp)
        tabl (i+1, 16) = HypCC(i).cell.(sprintf('Trace_%d',idx)).Freqs_after_Hyp(1);
    end
    
    % ---- ploting the trace ----- 
    subplot (nCells,2, nplot)
    hold on
        xlim ([0 10]);
        ylim([-0.1 0.05]);
        line ([0 10], [0 0], 'Color',[0.7 0.7 0.7]);
        plot (HypCC(i).cell.(sprintf('Trace_%d',idx)).trace(:,1), HypCC(i).cell.(sprintf('Trace_%d',idx)).trace(:,2));
        plot(HypCC(i).cell.(sprintf('Trace_%d',idx)).at52mv(1)/Fs,HypCC(i).cell.(sprintf('Trace_%d',idx)).at52mv(3),'Marker','x', 'Color', 'r');
        if isempty (HypCC(i).cell.(sprintf('Trace_%d',idx)).at45mv)
            plot(HypCC(i).cell.(sprintf('Trace_%d',idx)).at46mv(1)/Fs,HypCC(i).cell.(sprintf('Trace_%d',idx)).at46mv(3),'Marker','x', 'Color', 'r');
            line ([HypCC(i).cell.(sprintf('Trace_%d',idx)).at52mv(1)/Fs HypCC(i).cell.(sprintf('Trace_%d',idx)).at46mv(1)/Fs],...
                [HypCC(i).cell.(sprintf('Trace_%d',idx)).at52mv(3) HypCC(i).cell.(sprintf('Trace_%d',idxj)).at46mv(3)], 'Color', 'r');
        else
            plot(HypCC(i).cell.(sprintf('Trace_%d',idx)).at45mv(1)/Fs,HypCC(i).cell.(sprintf('Trace_%d',idx)).at45mv(3),'Marker','x', 'Color', 'r');
            line ([HypCC(i).cell.(sprintf('Trace_%d',idx)).at52mv(1)/Fs HypCC(i).cell.(sprintf('Trace_%d',idx)).at45mv(1)/Fs],...
                [HypCC(i).cell.(sprintf('Trace_%d',idx)).at52mv(3) HypCC(i).cell.(sprintf('Trace_%d',idx)).at45mv(3)], 'Color', 'r');
        end
        
        title ([HypCC(i).cell.name sprintf(' Trace_%d',idx)], 'Interpreter','none');
        set(gca,'FontSize', 4); 
    hold off
    nplot = nplot+1;
    
    subplot (nCells,2, nplot)
    ssb = 2;
    sse = 6;
    hold on
        xlim ([2 6]);
        ylim([-0.1 0.05]);
        line ([0 10], [0 0], 'Color',[0.7 0.7 0.7]);
        plot (HypCC(i).cell.(sprintf('Trace_%d',idx)).trace(ssb*Fs:sse*Fs,1), HypCC(i).cell.(sprintf('Trace_%d',idx)).trace(ssb*Fs:sse*Fs,2));
        plot(HypCC(i).cell.(sprintf('Trace_%d',idx)).at52mv(1)/Fs,HypCC(i).cell.(sprintf('Trace_%d',idx)).at52mv(3),'Marker','x', 'Color', 'r');
        if isempty (HypCC(i).cell.(sprintf('Trace_%d',idx)).at45mv)
            plot(HypCC(i).cell.(sprintf('Trace_%d',idx)).at46mv(1)/Fs,HypCC(i).cell.(sprintf('Trace_%d',idx)).at46mv(3),'Marker','x', 'Color', 'r');
            line ([HypCC(i).cell.(sprintf('Trace_%d',idx)).at52mv(1)/Fs HypCC(i).cell.(sprintf('Trace_%d',idx)).at46mv(1)/Fs],...
                [HypCC(i).cell.(sprintf('Trace_%d',idx)).at52mv(3) HypCC(i).cell.(sprintf('Trace_%d',idx)).at46mv(3)], 'Color', 'r');
        else
            plot(HypCC(i).cell.(sprintf('Trace_%d',idx)).at45mv(1)/Fs,HypCC(i).cell.(sprintf('Trace_%d',idx)).at45mv(3),'Marker','x', 'Color', 'r');
            line ([HypCC(i).cell.(sprintf('Trace_%d',idx)).at52mv(1)/Fs HypCC(i).cell.(sprintf('Trace_%d',idx)).at45mv(1)/Fs],...
                [HypCC(i).cell.(sprintf('Trace_%d',idx)).at52mv(3) HypCC(i).cell.(sprintf('Trace_%d',idx)).at45mv(3)], 'Color', 'r');
        end
    set(gca,'FontSize', 4);    
    hold off
    nplot=nplot+1;
    
end

filename = [GroupName '_CCHyp_closest-80mV.xlsx'];
xlswrite (filename, tabl);

fname = ['ALL_CChyp_raster_' GroupName];
saveas (f, fname, 'fig')
saveas (f, fname, 'svg')

%% --- plot ALL hyperpolarisation around -80 ON TOP of each other ---
if nCells > 1
    f = figure;
    set (f, 'Position', [50 50 1200 700]);
    meanTrace = [HypCC(1).cell.(sprintf('Trace_%d',1)).trace(2*Fs:6*Fs,1)];
    for i = 1:nCells % going through all sweeps/diferent stimulation
        j = HypCC(i).cell.nTraces;
        minV_SAGs = [];
        if j>1
            for k=1:j
                minV_SAGs = [minV_SAGs HypCC(i).cell.(sprintf('Trace_%d',k)).minSAG];
            end
            [n, idx] = min( abs(abs(minV_SAGs) - 80) );
        else idx = 1;
        end

        meanTrace = [meanTrace HypCC(i).cell.(sprintf('Trace_%d',idx)).trace(ssb*Fs:sse*Fs,2)];

        ssb = 2;
        sse = 6;
        hold on
            xlim ([2 6]);
            ylim([-0.1 0.05]);
            line ([0 10], [0 0], 'Color',[0.7 0.7 0.7]);
            plot (HypCC(i).cell.(sprintf('Trace_%d',idx)).trace(ssb*Fs:sse*Fs,1), ...
                HypCC(i).cell.(sprintf('Trace_%d',idx)).trace(ssb*Fs:sse*Fs,2), 'LineWidth', 0.05, 'Color',[0.3 0.3 0.3] );
            xlabel ('(s)');
            ylabel ('(mV)');

    end

    plot (meanTrace (:,1), mean(meanTrace(:,2:end)')', 'LineWidth', 3, 'Color',[0.1 0.1 0.1]);
     hold off
    fname = ['ALL_CChyp_overlapping_' GroupName];
    saveas (f, fname, 'fig')
    saveas (f, fname, 'svg')
end


%% export data tables as excel
% Cell names export
tabl4CellNames = repmat (blanks(12),nCells,1);
row = 1;
for i = 1:nCells
    NameLength = length (HypCC(i).cell.name);
    tabl4CellNames (i, 1:NameLength) = HypCC(i).cell.name; % 
end
     
fileID = fopen('CellName_closest-80mV.txt','wt');
for r=1:size(tabl4CellNames,1)
     fprintf(fileID,'%s\n',tabl4CellNames(r,:));
end
fclose(fileID);

disp('Done')

close all;