%% Clear the workspace

clear all
close all
clc

fclose('all');

%% Read in the 3T numbers and the BRUNS for the HCM genotype

[ Num, Txt, Raw ] = xlsread('GenScan II Subjects Directory with Genotype - Extended.xlsx', 'HCM');

Head = Raw(1, :);
Data = Raw(2:end, :);

ScannerID = Data(:, 1);
SubjectID = Data(:, 2);

ScannerID = cellfun(@strtrim, ScannerID, 'UniformOutput', false);
SubjectID = cellfun(@strtrim, SubjectID, 'UniformOutput', false);

Root = 'C:\Users\ptokarcz\Desktop\000 - Stable CSI Workflow\HCM';

Compilation = fullfile(Root, 'Results Compilation - HCM.pdf');

NFOLDERS = numel(ScannerID);

wb = waitbar(0, 'Creating HCM compilation');

for n = 1:NFOLDERS
  pft_CombinePdfOutputs(Root, ScannerID{n}, SubjectID{n}, Compilation);
  
  waitbar(double(n)/double(NFOLDERS), wb, sprintf('Processed %1d of %1d folders', n, NFOLDERS));
end

waitbar(1, wb, sprintf('Processed %1d of %1d folders', NFOLDERS, NFOLDERS));

pause(1.0);

delete(wb);

%% Read in the 3T numbers and the BRUNS for the TTNtv (High PSI) genotype

[ Num, Txt, Raw ] = xlsread('GenScan II Subjects Directory with Genotype - Extended.xlsx', 'TTNtv (High PSI)');

Head = Raw(1, :);
Data = Raw(2:end, :);

ScannerID = Data(:, 1);
SubjectID = Data(:, 2);

ScannerID = cellfun(@strtrim, ScannerID, 'UniformOutput', false);
SubjectID = cellfun(@strtrim, SubjectID, 'UniformOutput', false);

Root = 'C:\Users\ptokarcz\Desktop\000 - Stable CSI Workflow\TTNtv (High PSI)';

Compilation = fullfile(Root, 'Results Compilation - TTNtv (High PSI).pdf');

NFOLDERS = numel(ScannerID);

wb = waitbar(0, 'Creating TTNtv (High PSI) compilation');

for n = 1:NFOLDERS
  pft_CombinePdfOutputs(Root, ScannerID{n}, SubjectID{n}, Compilation);
  
  waitbar(double(n)/double(NFOLDERS), wb, sprintf('Processed %1d of %1d folders', n, NFOLDERS));
end

waitbar(1, wb, sprintf('Processed %1d of %1d folders', NFOLDERS, NFOLDERS));

pause(1.0);

delete(wb);

%% Read in the 3T numbers and the BRUNS for the TTNtv (PSI < 90) genotype

[ Num, Txt, Raw ] = xlsread('GenScan II Subjects Directory with Genotype - Extended.xlsx', 'TTNtv (PSI < 90)');

Head = Raw(1, :);
Data = Raw(2:end, :);

ScannerID = Data(:, 1);
SubjectID = Data(:, 2);

ScannerID = cellfun(@strtrim, ScannerID, 'UniformOutput', false);
SubjectID = cellfun(@strtrim, SubjectID, 'UniformOutput', false);

Root = 'C:\Users\ptokarcz\Desktop\000 - Stable CSI Workflow\TTNtv (PSI_LT_90)';

Compilation = fullfile(Root, 'Results Compilation - TTNtv (PSI_LT_90).pdf');

NFOLDERS = numel(ScannerID);

wb = waitbar(0, 'Creating TTNtv (PSI < 90) compilation');

for n = 1:NFOLDERS
  pft_CombinePdfOutputs(Root, ScannerID{n}, SubjectID{n}, Compilation);
  
  waitbar(double(n)/double(NFOLDERS), wb, sprintf('Processed %1d of %1d folders', n, NFOLDERS));
end

waitbar(1, wb, sprintf('Processed %1d of %1d folders', NFOLDERS, NFOLDERS));

pause(1.0);

delete(wb);

%% Read in the 3T numbers and the BRUNS for the Negative genotype

[ Num, Txt, Raw ] = xlsread('GenScan II Subjects Directory with Genotype - Extended.xlsx', 'Negative');

Head = Raw(1, :);
Data = Raw(2:end, :);

ScannerID = Data(:, 1);
SubjectID = Data(:, 2);

ScannerID = cellfun(@strtrim, ScannerID, 'UniformOutput', false);
SubjectID = cellfun(@strtrim, SubjectID, 'UniformOutput', false);

Root = 'C:\Users\ptokarcz\Desktop\000 - Stable CSI Workflow\Negative';

Compilation = fullfile(Root, 'Results Compilation - Negative.pdf');

NFOLDERS = numel(ScannerID);

wb = waitbar(0, 'Creating Negative compilation');

for n = 1:NFOLDERS
  pft_CombinePdfOutputs(Root, ScannerID{n}, SubjectID{n}, Compilation);
  
  waitbar(double(n)/double(NFOLDERS), wb, sprintf('Processed %1d of %1d folders', n, NFOLDERS));
end

waitbar(1, wb, sprintf('Processed %1d of %1d folders', NFOLDERS, NFOLDERS));

pause(1.0);

delete(wb);


  
  




