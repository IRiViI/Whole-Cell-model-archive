function bar_archive(archive, matrix, type, varargin)
% Present values in bar diagram form
%   bar_archive(archive, matrix, type)
%   archive: archive object
%   matrix: n x m matrix with all the values of a certain state file output
%   type. n: reactions/molecules, m number of data sets
%   type: Type of data (Metabolite, MetabolicReaction, ProteinComplex,...)
%
%   'select'     - Select which molecules should be included
%   'xAngle'     - Adjust the angle of the xLabel
%   'xLabel'     - The value of the info structure that should be displayed
%                  along the x axis.
%   'infoNumber' - Select info structure
%
%   Examples:
%       [y1,std1] = archive.average_value(1,100,'Metabolite','counts',1);
%       [y2,std2] = archive.average_value(2,100,'Metabolite','counts',1);
%       matrix = [y1,y2];
%       select = (abs(y2-y1)./std1>2);
%       archive.bar(matrix,'Metabolite','select',select);

% Default options
options.infoNumber = 1;             % The number of the info structure
options.angle = 45;                % Angle of label
options.xLabel = 'ID';              % Display the ID values along the x axis
options.select = 1:size(matrix,1);  % The molecules/reactions to include
options.preselected = false;        % If matrix is already trimmed down to the values which should be displayed
options.horizontal = false;

% Adjust options
inputOptions = struct(varargin{1:end});
fieldNames = fieldnames(inputOptions);
for i = 1:size(fieldNames,1)
    options.(fieldNames{i}) = inputOptions.(fieldNames{i});
end

% Convert select such that both type of inputs are valid. Example
% select = boolean([1,0,1]) or select = [1,3]
if islogical(options.select)
    options.select = find(options.select);
end

% Total number of elements selected
tSelect = length(options.select);

% values to present
if ~options.preselected
values = matrix(options.select,:);
else
    values = matrix;
end

% x values
x = 1:tSelect;

% Make bar diagram
label = {archive.info(options.infoNumber).(type)(options.select).(options.xLabel)};
if options.horizontal
barh(x, values, 'grouped')
% Adjust appearance
set(gca,'YTick',x);
set(gca,'YTickLabelRotation',options.angle);
set(gca,'YTickLabel',label);
else
bar(x, values, 'grouped')
% Adjust appearance
set(gca,'XTick',x);
set(gca,'XTickLabelRotation',options.angle);
set(gca,'XTickLabel',label);
end

end