%Extract data of a xls file and add it to the archive.
%
% archive = xls_to_archive(archive,'field1.field2.linker',...
%   'field1.field2.data','path/to/xls/file','xls_sheet_name',linker_column,
%   data_column);
%
%  In this file are the values in the column "data_column" of the
%  "xls_sheet_name" in the "path/to/xls/file" added to the "archive" in the
%  subfield "archive.field1.field2.data". The values are inserted in the order
%  according to the linkers.
%
% Example: element(1,42) is 'ATP' in the data sheet and
% archive.field1.field2.linker(33) or archive.field1.field2(33).linker is also
% 'ATP' then: archive.field1.field2.data(33) and archive.field1.field2(33).data
% are the value of element(4,42) respectively in the case where data_column
% is 4.
%
% Note:
% There are two different ways that the archive is structured. In one case is
% are all the values saved as one array under the subfield
% (archive.field.subfield = array). In the other case, are all the values
% saved under different entries of the field (archive.field(1).subfield =
% array(1), archive.field(2).subfield = array(2),...). Both cases are taken
% into account

% Author: Rick Vink, rickvink@mit.edu h.w.vink@student.tudelft.nl
% Affilitation: Timothy Lu, MIT
% Last updated: 04/12/2016

% archive=xls_to_archive(archive,'info.Metabolite.ID','info.Metabolite.Category',...
%             '/home/rick/Dropbox/Documents/Karr/mmc4.xlsx',...
%             'S3G-Metabolites',1,16);

function varargout = xls_to_archive(varargin)

% Mandatory inputs                  % Example:
archive = varargin{1};                 % archive
options.list_link = varargin{2};    % 'info.Metabolite.ID'
options.list_data = varargin{3};    % 'info.Metabolite.Category'
options.file_name = varargin{4};    % '/home/rick/Dropbox/Documents/Karr/mmc4.xlsx'
options.sheet_name = varargin{5};   % 'S3G-Metabolites'
options.file_link = varargin{6};    % 1
options.file_data = varargin{7};    % 16

% Default settings
options.field_name = options.sheet_name; % Default setting for field_name

% Change default settings
if nargin > 7
    input_options = struct(varargin(8:end));
    field_names = fieldnames(input_options);
    for i = 1:size(field_names,1)
        options.(field_names{i}) = input_options.(field_names{i});
    end
end

% ------------- Process inputs ------------- %

% ------------- list

% Get the fields of list to read and write
options.read_location = strsplit(options.list_link,'.');
options.write_location = strsplit(options.list_data,'.');

% Get the size of the read fields
structure = archive;
t_field = length(options.read_location);
for i_field = 1:t_field
    structure = structure.(options.read_location{i_field});
    options.read_length(i_field) = length(structure);
end

% ------------- Excel

% Get sheet names
[~,sheet_name] = xlsfinfo(options.file_name);

% Target sheet number
sheet_number = find(strcmp(sheet_name,options.sheet_name));

% Get the raw data
[~,~,sheet_raw] = xlsread(options.file_name,sheet_name{sheet_number});

% Number of rows
t_row = size(sheet_raw,1);

% ------------- Process requiest ------------- %

% If the values are individually saved under separate entries of the parent
% Then get the parent field
if options.read_length(end) < options.read_length(end-1)
    t_field = length(options.write_location);
    data = archive;
    for i_field = 1:t_field-1
        data = data.(options.write_location{i_field});
    end
elseif options.read_length(end) > options.read_length(end-1)
    data(1:options.read_length(end)) = {''};
end

% For every row in the data sheet, get data
for i_row = 1:t_row
    % Get the value and linker
    value = sheet_raw(i_row,options.file_data);
    link = sheet_raw(i_row,options.file_link);
    if options.read_length(end) > options.read_length(end-1)
        % If the values are saved as one array
        number = find(strcmp(archive.info.Metabolite.ID,link)); % Link the row with the archive
        if ~isempty(number)
            if isnan(value{1})
                data(number) = {''};
            else
                data(number) = value; % Set value
            end
        end
    elseif options.read_length(end) < options.read_length(end-1)
        % If the values are individually saved under separate entries of the parent
        number = find(strcmp({data.(options.read_location{end})},link)); % Link the row with the archive
        % Add value if match is found
        if ~isempty(number)
            data(number).(options.write_location{end}) = value{1}; % Set value
        end
    end
end

% Save Data

% Find the location that should be edited
if options.read_length(end) > options.read_length(end-1)
    replace_field = length(options.read_location);
elseif options.read_length(end) < options.read_length(end-1)
    replace_field = length(options.read_location) - 1;
end

% Save
a = options.write_location;
if replace_field == 0
    archive = data;
elseif replace_field == 1
    archive.(a{1}) = data;
elseif replace_field == 2
    archive.(a{1}).(a{2}) = data;
elseif replace_field == 3
    archive.(a{1}).(a{2}).(a{3}) = data;
elseif replace_field == 4
    archive.(a{1}).(a{2}).(a{3}).(a{4}) = data;
elseif replace_field == 5
    archive.(a{1}).(a{2}).(a{3}).(a{4}).(a{5}) = data;
else
    print('To large field structure, please edit function')
end

% Set output
varargout{1} = archive;

end