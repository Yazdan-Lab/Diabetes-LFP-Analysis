function  [SWR_files, LFP_files, SWRLTDIdx, rem] = initalize_animal(cur_animal, animal_list)

cd(animal_list(cur_animal).name)
% Get file names of specific animals relevant files
%%%%%%
load('SWR_Index.mat'); % load SWR HTD/LTD  indices
load('REM.mat'); % load HTD/LTD State times
SWR_files = dir('SWR_R_*'); % Grab the number of files that have SWR event timings (usually 2, sometimes 1 or 3)
SWR_files = {SWR_files.name}; % throw away useless info

LFP_files = dir('LFP*'); % grab number of LFP files (sometimes 1 or 3 as well)
LFP_files = {LFP_files.name}; % throw away useless info
end