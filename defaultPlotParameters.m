fontname = 'helvetica';
fontsize = 16;
ABCfontsize = 24;
fontweight = 'normal';
linewidth = 3;

matlabGrey = [1 1 1]*215/255;

global AZred AZblue AZcactus AZsky AZriver AZsand AZmesa AZbrick

AZred = [171,5,32]/256;
AZblue = [12,35,75]/256;
AZcactus = [92, 135, 39]/256;
AZsky = [132, 210, 226]/256;
AZriver = [7, 104, 115]/256;
AZsand = [241, 158, 31]/256;
AZmesa = [183, 85, 39]/256;
AZbrick = [74, 48, 39]/256;

% lighten blue
AZblue = 1*AZblue + 0*[1 1 1];

% rgb (171,5,32) AZ red
% rgb(12,35,75) blue
% rgb(92, 135, 39) cactus
% rgb(132, 210, 226) sky
% rgb(7, 104, 115) sand
% rgb(183, 85, 39) mesa
% rgb(74, 48, 39) brick
% 
% rgb(207, 224, 216) light cactus
% rgb(200, 217, 216) light sky
% rgb(182, 190, 193) light sand
% rgb(252, 225, 182) light mesa
% rgb(250, 231, 216) light brick
% rgb(230, 227, 217) light???


global orange
orange = [0.906 0.463 0.247];
% colormap gray
% CC = colormap;
% CM = (CC).*repmat((1-[0.906 0.463 0.247]/0.906), size(CC,1),1);
% colormap(1-CM)

set(0, 'defaultfigurecolor', 'w', ...
    'defaultaxesfontsize', fontsize, ...
    'defaultaxesfontweight', fontweight, ...
    'defaultaxesfontname', fontname, ...
    'defaultaxestickdir', 'out', ...
    'defaultaxesbox', 'off', ...
    'defaultaxesydir', 'normal', ...
    'defaultlinelinewidth', linewidth, ...
    'defaultlinemarkersize', 30, ...
    'defaultfigureposition', [811   486   618   500])

