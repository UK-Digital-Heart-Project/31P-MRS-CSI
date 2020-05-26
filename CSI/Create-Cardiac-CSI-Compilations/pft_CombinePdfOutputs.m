function pft_CombinePdfOutputs(Root, ScannerID, SubjectID, Compilation)

%% List the already existing PDF auditing files

A = fullfile(Root, ScannerID, 'Auditing Files', sprintf('%s - %s - CSI Matrix Plot.pdf', ScannerID, SubjectID));
B = fullfile(Root, ScannerID, 'Auditing Files', sprintf('%s - %s - Septal Spectra.pdf', ScannerID, SubjectID));
C = fullfile(Root, ScannerID, 'Auditing Files', sprintf('%s - %s - Averaged Septal Spectrum.pdf', ScannerID, SubjectID));

Inputs = { A, B, C };

%% Find the CSI Matrix Screenshots and convert them from PNG to PDF, however many there are

Listing = dir(fullfile(Root, ScannerID, '* CSI Matrix Screenshot *.png'));

Entries = { Listing.name };

NPIXELS = numel(Entries);

for n = 1:NPIXELS
  X = imread(fullfile(Root, ScannerID, Entries{n}));
  
  hf = figure('Name', 'CSI Matrix Screenshot', 'MenuBar', 'none', 'NumberTitle', 'off', 'Visible', 'off');
  ha = axes(hf);
  
  imshow(X);

  [ P, N, E ] = fileparts(Entries{n});
  
  Title = N;
  
  title(Title, 'Interpreter', 'none');

  pft_FormatLandscapePrinting(hf, ha);
  
  p = strfind(Entries{n}, '.');
  q = p - 1;
  F = sprintf('%s.pdf', Entries{n}(1:q));
  T = fullfile(Root, ScannerID, F);

  print(hf, T, '-dpdf', '-fillpage');
  
  Inputs = horzcat(Inputs, T);

  delete(ha);
  delete(hf);
end

%% Find the Amares Results Screenshot and convert that from PNG to PDF

S = fullfile(Root, ScannerID, sprintf('%s - %s - Amares Results Screenshot.png', ScannerID, SubjectID));
T = fullfile(Root, ScannerID, sprintf('%s - %s - Amares Results Screenshot.pdf', ScannerID, SubjectID));

Y = imread(S);

hf = figure('Name', 'Amares Results Screenshot', 'MenuBar', 'none', 'NumberTitle', 'off', 'Visible', 'off');
ha = axes(hf);

imshow(Y);

[ P, N, E ] = fileparts(S);
  
Title = N;
  
title(Title, 'Interpreter', 'none');

pft_FormatLandscapePrinting(hf, ha);
  
print(hf, T, '-dpdf', '-fillpage');

pause(0.25);
  
Inputs = horzcat(Inputs, T);

delete(ha);
delete(hf);

%% Combine the inputs to create a new and more comprehensive PDF summary file

Output = fullfile(Root, ScannerID, sprintf('%s - %s - Comprehensive Cardiac CSI Summary.pdf', ScannerID, SubjectID));

append_pdfs(Output, Inputs{:});

pause(0.25);

%% Update the PDF compilation

append_pdfs(Compilation, Output);

pause(0.25);

%% That's all !

end







