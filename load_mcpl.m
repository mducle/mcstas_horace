function out_struct = load_mcpl(filename)
% Function to load a Monte-Carlo Particle List binary file (MCPL)
% For more details see: https://mctools.github.io/mcpl/ 
% or https://doi.org/10.1016/j.cpc.2017.04.012
%
% Syntax:
% mcpl = load_mcpl(filename)
%
% The output structure has the following fields:
% .source    - A comment showing what generated this MCPL
% .comments  - A cell array of comments (converted to string)
% .blob_keys - A cell array of keys for user "blobs" (converted to string)
% .blobs     - A cell array of user "blobs" (also converted to string)
% .np        - Number of particles
% .pdg       - A scalar or np-element vector of the Monte Carlo Particle
%              Data Group (PDG) particle number scheme identifier
%              See: http://pdg.lbl.gov/mc_particle_id_contents.html
% .pos       - A 3 x np array of particle positions (in cm)
% .kin       - An np-element array of particle kinetic energy (in MeV)
% .dir       - A 3 x np array of (unit) particle momentum directions
% .time      - An np-element array of particle time (in milliseconds)
% .weight    - A scalar or np-element array of particle weights or intensity
%              If the data comes from a McStas simulation using the ISIS
%              moderator the intensity is in (n/nsim)/s/uA, where nsim
%              is the total number of neutrons used in the simulation.
% .pol       - (optional) A 3 x np array of particle polarisations.
% .user      - (optional) A np-element integer array of user-customised data

fid = fopen(filename);
  magi = char(fread(fid, 4, 'char')');
  if magi ~= 'MCPL';
      fclose(fid);
      error('Input file is not a MCPL');
  end
  ver = char(fread(fid, 3, 'char')');
  if ver ~= '003';
      fclose(fid);
      error(['Cannot read MCPL file format version "' ver '"']);
  end
  endian = lower(char(fread(fid, 1, 'char')));
  np = fread(fid, 1, 'uint64', endian);
  ncmts = fread(fid, 1, 'uint32', endian);
  nblobs = fread(fid, 1, 'uint32', endian);
  has_user = fread(fid, 1, 'uint32', endian);
  has_pol = fread(fid, 1, 'uint32', endian);
  is_single = fread(fid, 1, 'uint32', endian);
  if is_single
      data_precision = 'single';
      precision_bytes = 4;
  else
      data_precision = 'double';
      precision_bytes = 8;
  end
  pdg = fread(fid, 1, 'int32', endian);
  particle_data_length = fread(fid, 1, 'uint32', endian);
  has_universal_weight = fread(fid, 1, 'uint32', endian);
  if (has_universal_weight)
      outstruct.weight = fread(fid, 1, 'double', endian);
  end
  nal = fread(fid, 1, 'uint32', endian);
  % Explicitly assume comments are character arrays (strings)
  out_struct.source = char(fread(fid, nal, 'char', endian)');
  for icmts = 1:ncmts
      nal = fread(fid, 1, 'uint32', endian);
      out_struct.comments{icmts} = char(fread(fid, nal, 'char', endian)');
  end
  for iblob = 1:nblobs
      nal = fread(fid, 1, 'uint32', endian);
      % Explicitly assume keys are character arrays (strings)
      out_struct.blob_keys{iblob} = char(fread(fid, nal, 'char', endian)');
  end
  for iblob = 1:nblobs
      nal = fread(fid, 1, 'uint32', endian);
      % Explicitly assume blobs are character arrays (strings)
      out_struct.blobs{iblob} = char(fread(fid, nal, 'char', endian)');
  end
  % Work out the byte length of floating point data per particle
  nfields = 7;  % Always have 3 positions, 2 directions, energy and time
  if has_pol
      nfields = nfields + 3;
  end
  if ~has_universal_weight
      nfields = nfields + 1;
  end
  nskip = 0;
  if pdg == 0
      nskip = 4;   % Number of bytes to skip when reading floating data.
  end
  if has_user
      nskip = nskip + 4;
  end
  if (nfields * precision_bytes + nskip) ~= particle_data_length
      fclose(fid);
      error('Error parsing header - particle data lengths don''t match!');
  end
  out_struct.np = np;
  out_struct.pdg = pdg;
  % Now reads in the data-section in one go...
  data_position = ftell(fid);  % In case we have to rewind to read the ints (pdg/user)
  data = fread(fid, np * nfields, data_precision, nskip, endian);  % Floating point data
  pos_ind = 1;  % Matlab, like Fortran, index from 1
  if has_pol
      out_struct.pol = [data(1:nfields:end) data(2:nfields:end) data(3:nfields:end)];
      pos_ind = 4;
  end
  out_struct.pos = [data((pos_ind):nfields:end) ...
                    data((pos_ind+1):nfields:end) ...
                    data((pos_ind+2):nfields:end)];
  fp1 = data((pos_ind+3):nfields:end);
  fp2 = data((pos_ind+4):nfields:end);
  out_struct.kin = data((pos_ind+5):nfields:end);
  % Do the inverse transform
  sig1 = fp1 <= 1;  % The <= (\leq) is needed because we use the negation
  sig2 = fp2 <= 1;  %   operator (~) later to get fp1 > 1 or fp2 > 1.
  % Note that Matlab does not allow user code to handle -0 (it maps it to +0)
  % The sign function below returns 1 for >0, -1 for <0 and 0 for -0 or +0
  signu = sign(out_struct.kin); 
  signu(signu==0) = 1;    % Convert 0's to 1's (assume positive)
  ux = zeros(size(fp1));
  uy = zeros(size(fp1));
  uz = zeros(size(fp1));
  case1 = find(~sig1 .* sig2);
  case2 = find(sig1 .* ~sig2);
  case3 = find(sig1 .* sig2);
  ux([case2; case3]) = fp1([case2; case3]);
  uy([case1; case3]) = fp2([case1; case3]);
  uz(case1) = 1 ./ fp1(case1);
  uz(case2) = 1 ./ fp2(case2);
  ux(case1) = signu(case1) .* sqrt(1 - uy(case1).^2 - uz(case1).^2);
  uy(case2) = signu(case2) .* sqrt(1 - ux(case2).^2 - uz(case2).^2);
  uz(case3) = signu(case3) .* sqrt(1 - ux(case3).^2 - uy(case3).^2);
  out_struct.dir = [ux uy uz];
  out_struct.time = data((pos_ind+6):nfields:end);
  if ~has_universal_weight
      out_struct.weight = data((pos_ind+7):nfields:end);
  end
  if nskip > 0
      fseek(fid, data_position + (nfields * precision_bytes), 'bof');
      if pdg == 0
          out_struct.pdg = fread(fid, np, 'int32', (nfields * precision_bytes + nskip - 4), endian);
          fseek(fid, data_position + (nfields * precision_bytes), 'bof');
      end
      if has_user
          out_struct.user = fread(fid, np, 'uint32', (nfields * precision_bytes + nskip - 4), endian);
      end
  end
fclose(fid);