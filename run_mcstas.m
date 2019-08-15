function out = run_mcstas(instrument, ei, freq, chopper, nsim, output_file, sample)
% Runs a McStas simulation for a given Ei, freq, chopper and returns the
% resulting neutron trajectory list at the sample position.
% 
% out = run_mcstas(instrument, ei, freq, chopper, nneut, output_file, sample)
%
% instrument - instrument name. 'map[s]', 'mar[i]', 'me[rlin]', 'l[et]'
%              are recognised (only part before brackets are needed)
% ei - incident energy (meV)
% freq - frequency of the Fermi or resolution choppers (Hz)
% chopper - chopper type ('a', 'b', 'c', 'r', 's', 'g', 'high flux',
%           'high resolution', and 'intermediate' are recognised).
% nsim - number of neutron trajectories to generate at moderator. 
%        (around an order of magnitude of 1% get through to sample).
%        (default: 1e7 ~ 2min on a modern i7 [single-core])
% output_file - the output file name (default: 'mcstas.mcpl'). If this file
%               exists, the function will try to check if the ei, freq, etc
%               matches, and if so will load from the existing file.
% sample - either "cylinder" (default for MARI) or "plate" (for others)
%
% out - is a struct which takes the format of the load_mcpl function.
%       in particular, the data fields are:
%       .pos - the (x,y,z) position in cm.
%       .kin - the kinetic energy in MeV (mega-electron-volts)
%       .dir - the fractional directional (unit velocity) vector
%       .time - the ToF of the neutron in ms (miliseconds)
%       .weight - the weight of this trajectory in (n/nsim)/s/uA

if ~exist('output_file', 'var')
    output_file = 'mcstas.mcpl';
end

exec_path = strrep(mfilename('fullpath'), mfilename, '');
setenv('MCTABLES', exec_path);

inst = lower(instrument);
if strncmp(inst, 'map', 3)
    mcstas_exec = 'MAPS_GuideAsTender';
    def_samp = 'plate';
elseif strncmp(inst, 'mar', 3)
    mcstas_exec = 'ISIS_MARI_upgraded';
    def_samp = 'cylinder';
elseif strncmp(inst, 'me', 2)
    mcstas_exec = 'merlin_RAEmod';
    def_samp = 'plate';
elseif strncmp(inst, 'l', 1)
    mcstas_exec = 'let_new';
    def_samp = 'plate';
else
    error(sprintf('Unrecognised instrument name: %s', instrument));
end

if ~exist('sample', 'var')
    sample = def_samp;
end

if exist(output_file, 'file')
    out = load_mcpl(output_file);
    idx = find(cellfun(@(x)strcmp(x, 'mccode_cmd_line'), out.blob_keys));
    if ~isempty(idx)
        rem = out.blobs{idx}; 
        ii=1; 
        while length(rem)>0
            [str{ii}, rem] = strtok(rem, '= '); 
            ii=ii+1; 
        end
        inst_saved = str{1};
        ei_saved = str2num(str{find(cellfun(@(x)strcmp(x, 'Ei'), str))+1});
        freq_saved = str2num(str{find(cellfun(@(x)strcmp(x, 'freq'), str))+1});
        chopper_saved = str{find(cellfun(@(x)strcmp(x, 'chopper'), str))+1};
        sample_saved = str{find(cellfun(@(x)strcmp(x, 'sample'), str))+1};
        if contains(inst_saved, mcstas_exec) && ei_saved==ei && freq_saved==freq ...
            && strcmp(chopper_saved, chopper) && strcmp(sample_saved, sample)
            return
        end
    end
end

cmdstr = sprintf('%s\\%s -n%d Ei=%f freq=%d chopper="%s" output_filename=%s sample=%s', ...
    exec_path, mcstas_exec, nsim, ei, freq, chopper, output_file, sample);
[status, cmdout] = system(cmdstr);
if status == 0
    out = load_mcpl(output_file);
else
    error(sprintf('McStas failed, with this message:\n%s', cmdout));
end