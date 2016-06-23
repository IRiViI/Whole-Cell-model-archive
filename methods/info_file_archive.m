function info = info_file_archive(archive, iSet, iSimulation, file, varargin)
% Get information about file (date, storage space, directory status,
% datenum)
%   info = info_file_archive(archive, set, simulation, file)

% Folder of simulation
folder = archive.set(iSet).simulation(iSimulation).folder;

% Get file info of file
info = dir([folder '/' file]);

end