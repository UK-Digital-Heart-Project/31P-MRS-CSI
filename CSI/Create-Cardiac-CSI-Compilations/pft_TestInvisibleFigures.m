%% Clear the workspace

clear all
close all
clc

fclose('all');

% List the figures to be concatenated

Listing = dir(fullfile(pwd, '*Amares Results Screenshot.png'));

Entries = { Listing.name };

x = imread(Entries{1});

hf = figure('Name', 'Amares Results Screenshot', 'MenuBar', 'none', 'NumberTitle', 'off', 'Visible', 'off');
ha = axes(hf);

imshow(x);

title('3T7109 - 14TD01766 - Amares Results Screenshot');

pft_FormatLandscapePrinting(hf, ha);

print(hf, 'Amares Results Screenshot.pdf', '-dpdf', '-fillpage');

delete(ha);
delete(hf);





