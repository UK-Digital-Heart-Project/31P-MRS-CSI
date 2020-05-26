%% Clear the workspace

clear all
close all
clc

fclose('all');

%% Read in the 3T numbers and the BRUNS

% [ Num, Txt, Raw ] = xlsread('GenScan II Subjects Directory with Genotype.xlsx', 'GenScan II Subjects');
% 
% Head = Raw(1, :);
% Data = Raw(2:end, :);
% 
% ScannerID = Data(:, 1);
% SubjectID = Data(:, 2);
% 
% ScannerID = cellfun(@strtrim, ScannerID, 'UniformOutput', false);
% SubjectID = cellfun(@strtrim, SubjectID, 'UniformOutput', false);

ScannerID = { '3T7463';    '3T7466';    '3T7470' };
SubjectID = { '14JE02931'; '14JK01282'; '14JH01338' };
    
% Root = 'C:\Users\ptokarcz\Desktop\31-P March 2020\Early 2020\Source Data\Cardiac CSI';
  Root = 'C:\Users\ptokarcz\Desktop\TestTestTest';

Compilation = fullfile(Root, 'Cardiac CSI Results Compilation - Early 2020.pdf');

NFOLDERS = numel(ScannerID);

wb = waitbar(0, 'Creating cardiac CSI compilation');

for n = 1:NFOLDERS
  pft_CombinePdfOutputs(Root, ScannerID{n}, SubjectID{n}, Compilation);
  
  waitbar(double(n)/double(NFOLDERS), wb, sprintf('Processed %1d of %1d folders', n, NFOLDERS));
end

waitbar(1, wb, sprintf('Processed %1d of %1d folders', NFOLDERS, NFOLDERS));

pause(1.0);

delete(wb);