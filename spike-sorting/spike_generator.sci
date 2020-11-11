// Created by Robert Toth, July 2018

function spike_generator(path)

  // Sanitize input path
  path=fullpath(path);

  // Read user params
  params=input_parser(path);

  // Generate spike times
  pos=timepoints(params);

  // Read and normalize template file if user params require
  templates=read_templates(path,params);

  // Display progress bar
  winH=waitbar('Generating spike data');

  // Case 1 - Amplitude Sweep
  if params.mode==1 then

    // Find number of items to generate for progress bar
    waitmax=(params.amplitudeSweep1*params.amplitudeSweep2*params.repetitions);

    for i=1:params.repetitions
      noise=grand(params.interval,params.channels,'nor',0,params.noiseAmplitude);

      counter=0;

      for j=1/params.amplitudeSweep1:1/params.amplitudeSweep1:1
        for k=1/params.amplitudeSweep2:1/params.amplitudeSweep2:1
          // Scale template amplitudes
          [temp,a1,a2]=scale_templates(templates,params,j,k);

          // Add spikes at the generated positions
          spikes=zeros(params.interval,params.channels);
          for l=1:params.numberOfSpikes
            clu(l,1)=modulo(l,params.spikeTypes);
            spikes(pos(l):(pos(l)+params.spikeLength-1),:) = ...
            temp(:,(clu(l,1)*params.channels+1):((clu(l,1)+1)*params.channels));
          end

          signal=noise+spikes;

          counter=counter+1;

          spike_export(path,['Set_'+string(i),string(counter)],params,signal,pos,clu+1);

          amp1(1,counter)=a1;
          amp2(1,counter)=a2;

          waitbar(((i-1)*params.amplitudeSweep1*params.amplitudeSweep2+counter)/waitmax,winH);
        end
      end
    end
    // Save the amplitudes
    csvWrite(matrix(amp1,params.amplitudeSweep1,params.amplitudeSweep2),fullfile(path,'amp1.csv'),',','.');
    csvWrite(matrix(amp2,params.amplitudeSweep1,params.amplitudeSweep2),fullfile(path,'amp2.csv'),',','.');

    // Case 2 - Channel Combinations
  elseif params.mode==2 then

    // Find number of items to generate for progress bar
    waitmax=0;
    for i=1:params.channels
      waitmax=waitmax+size(subsets(params.channels,i),1);
    end
    waitmax=waitmax*params.repetitions;

    counter=0;

    probe=[1,8,2,7,3,6,4,5]; // Lengthwise order of probe sites on Buzsaki32

    for i=1:params.repetitions
      for j=1:params.channels
        temp=template_generator(templates,params,j);
        params.channelsInUse=probe((params.channels-j+1):1:$);
        noise=grand(params.interval,j,'nor',0,params.noiseAmplitude);

        for k=1:size(temp,3)
          // Add spikes at the generated positions
          spikes=zeros(params.interval,j);
          for l=1:params.numberOfSpikes
            clu(l,1)=modulo(l,params.spikeTypes);
            spikes(pos(l):(pos(l)+params.spikeLength-1),:) = temp(:,(clu(l,1)*j+1):((clu(l,1)+1)*j),k);
          end

          signal=noise+spikes;

          spike_export(path,['Set_'+string(i),string(j)+'_channels',string(k)],params,signal,pos,clu+1);

          counter=counter+1;
          waitbar(counter/waitmax,winH);
        end
      end
    end
  end

  close(winH);

endfunction


// Reads user parameters from params.txt into a struct
function params=input_parser(path)

  try
    // fullfile() constructs a file path appropriate for the user's operating system
    filename=fullfile(path,'params.txt');
    temp=csvRead(filename,'=',[],'string');
    var=temp(:,1); // Variable names
    val=evstr(temp(:,2)); // Values

    params.channels=val(find(var=='channels')); // Number of channels for template file
    params.recordingRate=val(find(var=='recordingRate')); // Sampling rate in [Hz]
    params.recordingLength=val(find(var=='recordingLength')); // Length of recording in [ms]
    params.spikeTypes=val(find(var=='spikeTypes')); // Number of spike types in template file
    params.spikeLength=val(find(var=='spikeLength')); // Length of a spike in [samples]
    params.spikesPerType=val(find(var=='spikesPerType')); // Number of spikes to generate for each type
    params.maxAmplitude=val(find(var=='maxAmplitude')); // If not NaN, spikes primary amplitudes normalized
    params.noiseAmplitude=val(find(var=='noiseAmplitude')); // Standard deviation of Gaussian noise
    params.repetitions=val(find(var=='repetitions')); // Generate multiple datasets
    params.completeness=val(find(var=='completeness')); // Required cluster completeness
    params.purity=val(find(var=='purity')); // Required cluster purity
    params.tolerance=val(find(var=='tolerance')); // Temporal tolerance for spike matching [samples]
    params.penalty=val(find(var=='penalty')); // Parameter k, 0- penalty for spike sorting errors 1- no penalty
    params.mode=val(find(var=='mode')); // Data generation mode. 1- Amplitude sweep, 2- Channel permutations

    // Generate derivative parameters
    params.minDistance=2*params.spikeLength; // Minimum samples between two spikes
    params.interval=params.recordingLength*params.recordingRate; // Length of recording in [samples]
    params.numberOfSpikes=params.spikeTypes*params.spikesPerType; // Total number of spikes

    if params.mode==1 then
      params.channelsInUse=1:1:params.channels;
      params.amplitudeSweep1=val(find(var=='amplitudeSweep1')); // Points for primary sweep
      params.amplitudeSweep2=val(find(var=='amplitudeSweep2')); // Points for secondary sweep
    end

  catch
    error('Parameters file malformatted or does not exist');

  end

endfunction


//Generates timepoints without overlap in 'interval'
function pos=timepoints(params)

  test = round((params.interval-(params.spikeLength-1))*grand(1, params.numberOfSpikes*5, 'def'));
  pos = test(1);
  counter = 2;

  for j = 2:size(test,2)
    // Choose a point to test
    current = test(j);
    // Check distance from stored points
    distances = abs(current - pos);
    if min(distances) >= params.minDistance then
      pos(counter) = current;
      counter = counter + 1;
    end
    if counter == params.numberOfSpikes+1 then
      break;
    end
  end

  if size(pos,1)<params.numberOfSpikes then
    disp('Recording length too short for random point generation, distributing spike times uniformly');
    pos=[params.spikeLength:...
    floor((params.interval-2*params.spikeLength)/params.numberOfSpikes):...
    params.interval-2*params.spikeLength]';
    pos=pos(grand(1, 'prm', (1:params.numberOfSpikes))); // Points are mixed to randomize type assignments
  end

endfunction


// Read and normalize template file if user params require
function templates=read_templates(path,params)

  // Read set of test spikes
  // fullfile() constructs a file path appropriate for the user's operating system
  filename=fullfile(path,'templates.csv');
  templates=csvRead(filename,',','.','double');

  // For handling positive and negative spikes
  //if max(max(templates,'r'),'c')>max(max(-1*templates,'r'),'c') then
  //    templates=-1*templates;
  //end

  // If params.maxAmplitude is defined, scale the templates
  if ~isnan(params.maxAmplitude) then
    for i=0:params.spikeTypes-1
      temp=templates(:,params.channels*i+1:params.channels*(i+1));
      scale=abs(params.maxAmplitude/min(min(temp,'r'),'c'));
      templates(:,params.channels*i+1:params.channels*(i+1))=temp.*scale;
    end
  end

endfunction


// Scales spike templates for amplitude sweep simulations
function [newtemplates,amp1,amp2]=scale_templates(templates,params,f1,f2)

  newtemplates=zeros(params.spikeLength,params.channels*params.spikeTypes);

  temp=zeros(params.spikeLength,params.channels);

  for i=1:params.spikeTypes
    temp=templates(:,((i-1)*params.channels+1):(i*params.channels));
    [%,ord]=gsort(min(temp,'r'),'g','i'); // Get amplitude order of the spike over all channels
    s=f2*f1*min(temp(:,ord(1,1)))/min(temp(:,ord(1,2))); // Scaling of secondary channels
    temp(:,ord(1,1))=f1*temp(:,ord(1,1)); // Scale the largest channel
    temp(:,setdiff(1:params.channels,ord(1,1)))=s*temp(:,setdiff(1:params.channels,ord(1,1))); // Scale all other channels
    newtemplates(:,((i-1)*params.channels+1):(i*params.channels))=temp;
    a1(i,1)=-min(temp(:,ord(1,1)));
    a2(i,1)=-min(temp(:,ord(1,2)));
  end

  amp1=mean(a1,1);
  amp2=mean(a2,1);

endfunction


// Exports spike data for KlustaKwik input format. Generates KlustaKwik settings
// files using a reduced Buzsaki32 probe. Exports spike classes and spike times
function spike_export(path,newfolders,params,signal,pos,clu)

  for i=1:size(newfolders,2)
    mkdir(path,newfolders(1,i));
    path=fullfile(path,newfolders(1,i));
  end

  [pos,order]=gsort(pos,'r','i');

  // fullfile() constructs a file path appropriate for the user's operating system
  posfilename=fullfile(path,'pos.csv');
  clufilename=fullfile(path,'clu.csv');

  csvWrite(pos,posfilename,',','.');
  csvWrite(clu(order,1),clufilename,',','.');

  prm_generator(path,params);
  prb_generator(path,params);

  // fullfile() constructs a file path appropriate for the user's operating system
  h5filename=fullfile(path,'test_signal.raw.kwd');
  h5file = h5open(h5filename, 'a');
  h5group(h5file, '/recordings');
  h5group(h5file, '/recordings/0')
  h5write(h5file, '/recordings/0/data', int16(signal'), 'H5T_STD_I16LE');
  h5close(h5file);

endfunction


// Generates a klustakwik .prm file
function prm_generator(path,params)

  // Defines content of prm file, ascii(10) used for line breaks
  prm=strcat([...
  '# pay attention to use '' and not ‘',ascii(10),...
  ascii(10),...
  'experiment_name = ''test_signal''',ascii(10),...
  'prb_file = ''probe.prb''',ascii(10),...
  ascii(10),...
  'traces = dict(',ascii(10),...
  '    raw_data_files=''test_signal.raw.kwd'',',ascii(10),...
  '    voltage_gain=10.,',ascii(10),...
  '    sample_rate=',string(params.recordingRate),',',ascii(10),...
  '    n_channels=',string(length(params.channelsInUse)),',',ascii(10),...
  '    dtype=''int16'',',ascii(10),...
  ')',ascii(10),...
  ascii(10),...
  'spikedetekt = dict(',ascii(10),...
  '    filter_low=500.,  # Low pass frequency (Hz)',ascii(10),...
  '    filter_high_factor=0.95 * .5,',ascii(10),...
  '    filter_butter_order=3,  # Order of Butterworth filter.',ascii(10),...
  ascii(10),...
  '    filter_lfp_low=0,  # LFP filter low-pass frequency',ascii(10),...
  '    filter_lfp_high=300,  # LFP filter high-pass frequency',ascii(10),...
  ascii(10),...
  '    chunk_size_seconds=1,',ascii(10),...
  '    chunk_overlap_seconds=.015,',ascii(10),...
  ascii(10),...
  '    n_excerpts=50,',ascii(10),...
  '    excerpt_size_seconds=1,',ascii(10),...
  '    threshold_strong_std_factor=4.5,',ascii(10),...
  '    threshold_weak_std_factor=2.,',ascii(10),...
  '    detect_spikes=''negative'',',ascii(10),...
  ascii(10),...
  '    connected_component_join_size=1,',ascii(10),...
  ascii(10),...
  '    extract_s_before=',string(floor(params.spikeLength/2)),',',ascii(10),...
  '    extract_s_after=',string(floor(params.spikeLength/2)),',',ascii(10),...
  ascii(10),...
  '    n_features_per_channel=3,  # Number of features per channel.',ascii(10),...
  '    pca_n_waveforms_max=10000,',ascii(10),...
  ')',ascii(10),...
  ascii(10),...
  'klustakwik2 = dict(',ascii(10),...
  '    num_starting_clusters=100,',ascii(10),...
  ')',ascii(10)...
  ]);

  // Create output file, write data
  // fullfile() constructs a file path appropriate for the user's operating system
  filename=fullfile(path,'params.prm');
  fd = mopen(filename,'wt');
  mfprintf(fd,prm);
  mclose(fd);

endfunction


function prb_generator(path,params)

  // Adjacency matrix of desired probe, here a single shank of Buzsaki32
  // Note, that Scilab uses 1-based indexing, while Klusta uses 0-based
  adj=zeros(8,8);
  adj(1,2)=1; adj(1,8)=1;
  adj(2,3)=1; adj(2,7)=1; adj(2,8)=1;
  adj(3,4)=1; adj(3,6)=1; adj(3,7)=1;
  adj(4,5)=1; adj(4,6)=1;
  adj(5,6)=1;
  adj(6,7)=1;
  adj(7,8)=1;

  adj=adj';

  // Remove unused channels from adjacency matrix
  if params.mode==2 then
    if length(params.channelsInUse)==1 then
      adj=1;
    else
      rm=setdiff(1:params.channels,params.channelsInUse); // The unused channels
      if length(rm)~=0 then
        adj(rm,:)=[];
        adj(:,rm)=[];
      end
    end
  end

  // Create a formatted list of adjacency edges
  adjacencyList='';
  for i=1:length(params.channelsInUse)
    for j=1:i
      if adj(i,j)==1 then
        adjacencyList=strcat([adjacencyList,ascii(10),...
        '                (',string(j-1),',',string(i-1),'),']); // -1 for indexing mismatch
      end
    end
  end

  // Construct the probe file
  prb=strcat([...
  'channel_groups = {',ascii(10),...
  '    # Shank index.',ascii(10),...
  '    0:',ascii(10),...
  '        {',ascii(10),...
  '            # List of channels to keep for spike detection.',ascii(10),...
  '            ''channels'': range(',string(length(params.channelsInUse)),'),',ascii(10),...
  '            ',ascii(10),...
  '            # Adjacency graph. Dead channels will be automatically discarded',ascii(10),...
  '            # by considering the corresponding subgraph.',ascii(10),...
  '            ''graph'': [',...
  adjacencyList,ascii(10),...
  '            ],',ascii(10),...
  '    }',ascii(10),...
  '}',ascii(10)...
  ]);

  // Create output file, write data
  // fullfile() constructs a file path appropriate for the user's operating system
  filename=fullfile(path,'probe.prb');
  fd = mopen(filename,'wt');
  mfprintf(fd,prb);
  mclose(fd);

endfunction

// Generates k element permutations of a template file
function out=template_generator(source,params,k)

  y=subset(params.channels,k);

  out=zeros(params.spikeLength,params.channels*k,size(y,1));

  for i=1:size(y,1)
    for j=1:k
      out(:,j:k:((size(source,2)/params.channels)*k),i)=source(:,y(i,j):params.channels:$);
    end
  end

endfunction


//The function generates m element subsets of set 1:n
//Source: https://groups.google.com/forum/#!msg/comp.soft-sys.math.scilab/HyYwDyGeT5U/A_0u4YjDNlEJ
//Available in Matlab as nchoosek(1:n,m)
//Available in Scilab 6.0.x only with atomsInstall("specfun"), specfun_subset(1:n,m)
function [A] = subset(n,m)
  ncT = factorial(n)/(factorial(m)*factorial(n-m));
  A(ncT,m) = 0;
  for j = 1:m
    k=j;
    i0=1;
    while i0<=ncT
      n0 = n-k; m0 = m-j;
      i1 = i0-1 + factorial(n0)/(factorial(m0)*factorial(n0-m0));
      A(i0:i1,j)=k;
      if(n-k == m-j & i0<ncT) then
        k = A(i1+1,j-1)+1;
      else
        k=k+1;
      end
      i0=i1+1;
    end
  end
endfunction

