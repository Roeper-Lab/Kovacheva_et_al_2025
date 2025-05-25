%% Spiketrain_Import
% By Kau? Machado Costa

% This script imports spike time data (spike times extracted from
% the Igor analysis) from .txt files

%%%Ideas:
% 
%

% Last modified on 25.04.2016

%% How to load original data
% You should load a .txt file with the spike times. Each neuron's data
% should be on its own column and the first line represents the column
% header (colheaders; ideally the documentation code for each neuron).
% Columns should be separated by tabs and empty spaces should NOT be filled
% with any data (Matlab will automatically fill them with "NAN"s).

%% Parameter definition, data import and structure formatting

clear all
GroupName = 'OHDA_late'; % ACSFmed_early

scale = 's';

% if strcmp(scale, 's') == 1;
%     binsize = 0.01;
%     xvalues = 0:binsize:45;
%     xadd = 0.2;
%     refpd = 0.005;
% elseif strcmp(scale, 'ms') == 1;
%     xvalues = 0:5:20^3;
%     xadd = 200;
%     refpd = 2.5;
% end


import = importdata('in_vivo_late_6OHDA_ISIs.txt');
colheaders = import.colheaders;
TimeStamps = import.data;
clear import

CellLabels = colheaders;

nCells = length(CellLabels);                                                % Gets number of cells
for c = 1:nCells
    Data(c).(GroupName).CellName = CellLabels (:,c);
    
    temp = TimeStamps(:,c);
    temp = temp(~isnan(temp));
    
    Data(c).(GroupName).TimeStamps = temp;
    
    Data(c).(GroupName).ISI = diff(Data(c).(GroupName).TimeStamps);         % Extracts ISIs
    Data(c).(GroupName).Period = 'cell';
    clear temp
end

% Set folders
mkdir(['ISI_' char(GroupName)]);
oldFolder = cd(['ISI_' char(GroupName)]);      % Sets the return to old directory

close all
clearvars -except Data nCells GroupName oldFolder;
save (['ISI_' char(GroupName)]);
cd(oldFolder)
clearvars