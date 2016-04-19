% Compress output values of the simulations
%
%  save_figure(figNumber,outputDir,fileName,type)
%  Example: save_figure(1,'path/to/output/folder',{'pdf','fig'})
%  save_figure(figNumber,outputDir,fileName)
%  filetype is pdf and fig (type = {'pdf','fig'}).

% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT, TU Deflt
% Last updated: 17/04/2016

function save_figure_archive(varargin)


% Process inputs
archive = varargin{1};
options.name = varargin{2};

options.dir = [archive.settings.dir '/' 'output'];
options.fit = 'On';
options.type = {'pdf','fig','png'};

for i = 1:length(varargin)
    if strcmp(varargin{i},'type')
      options.type = varargin{i+1}; 
   end
   if strcmp(varargin{i},'fit')
      options.fit = varargin{i+1}; 
   end
   if strcmp(varargin{i},'dir')
      options.dir = varargin{i+1}; 
   end
end

% Adjust figure 
if strcmp(options.fit,'On')
    set(gca,'units','centimeters')
    pos = get(gca,'outerPosition');
    set(gcf, 'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [pos(3) pos(4)]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition',[0 0 pos(3) pos(4)]);
else
    set(gcf, 'PaperPositionMode','auto');
end

% Make folder
[~,~,~] = mkdir(options.dir);

% Save files
for i = 1:size(options.type,2)                  % Save for every file type
    fullOutput = [options.dir '/' options.name '.' options.type{i}];
    saveas(gcf,fullOutput);                       % Save figure
% savefig(fullOutput);
end

end
