// Created by Robert Toth, July 2018

function sorting_quality(path)

  // Sanitize input path
  path=fullpath(path);

  // Read user params
  params=input_parser(path);

  winH=waitbar('Processing Klustakwik results');

  // Case 1 - Amplitude Sweep
  if params.mode==1 then

    // Find number of items to process for progress bar
    waitmax=(params.amplitudeSweep1*params.amplitudeSweep2*params.repetitions);

    count=0;

    for i=1:params.repetitions
      for j=1:params.amplitudeSweep1*params.amplitudeSweep2

        pfile=fullfile(path,'Set_'+string(i),string(j),'pos.csv');
        gpos=csvRead(pfile,',',[],'double');

        cfile=fullfile(path,'Set_'+string(i),string(j),'clu.csv');
        gclu=csvRead(cfile,',',[],'double');

        kwikfile=fullfile(path,'Set_'+string(i),string(j),'test_signal.kwik');
        kwik=h5open(kwikfile,'r');
        apos=double(h5read(kwik,'/channel_groups/0/spikes/time_samples'))';
        aclu=double(h5read(kwik,'/channel_groups/0/spikes/clusters/main'))';
        h5close(kwik);

        conf=confusion(params,gpos,gclu,apos,aclu);

        n(i,j)=sum((completeness(conf)>=params.completeness)&(purity(conf)>=params.purity));

        // DEBUG USE: output the confusion matrices
        // csvWrite(conf,fullfile(path,'Set_'+string(i),'conf'+string(j)+'.csv'),',','.');

        count=count+1;
        waitbar(count/waitmax,winH);
      end
    end
    cpc=matrix(mean(n,1),params.amplitudeSweep1,params.amplitudeSweep2)./params.channels;
    csvWrite(cpc,fullfile(path,'cpc.csv'),',','.');

    // Case 2 - Channel Combinations
  elseif params.mode==2 then

    // Find number of items to process for progress bar
    waitmax=0;
    for i=1:params.channels
      waitmax=waitmax+size(subsets(params.channels,i),1);
    end
    waitmax=waitmax*params.repetitions;

    n=zeros(params.repetitions,params.channels);
    count=0;

    for i=1:params.repetitions
      for j=1:params.channels
        for k=1:size(subsets(params.channels,j),1)
          pfile=fullfile(path,'Set_'+string(i),string(j)+'_channels',string(k),'pos.csv');
          gpos=csvRead(pfile,',',[],'double');

          cfile=fullfile(path,'Set_'+string(i),string(j)+'_channels',string(k),'clu.csv');
          gclu=csvRead(cfile,',',[],'double');

          kwikfile=fullfile(path,'Set_'+string(i),string(j)+'_channels',string(k),'test_signal.kwik');
          kwik=h5open(kwikfile,'r');
          apos=double(h5read(kwik,'/channel_groups/0/spikes/time_samples'))';
          aclu=double(h5read(kwik,'/channel_groups/0/spikes/clusters/main'))';
          h5close(kwik);

          conf=confusion(params,gpos,gclu,apos,aclu);

          // DEBUG USE: output the confusion matrices
          // csvWrite(conf,fullfile(path,'Set_'+string(i),string(j)+'_channels','conf'+string(k)+'.csv'),',','.');

          n0=sum((completeness(conf)>=params.completeness)&(purity(conf)>=params.purity));
          n(i,j)=n(i,j)+(n0./size(subsets(params.channels,j),1))./j;

          count=count+1;
          disp(count);
          waitbar(count/waitmax,winH);
        end
      end
    end

    // DEBUG USE: output cluster per channel for each trial
    csvWrite(n,fullfile(path,'cpc_full.csv'),',','.');

    cpc(1,:)=mean(n,1);
    cpc(2,:)=stdev(n,1);
    csvWrite(cpc,fullfile(path,'cpc.csv'),',','.');

  else
    error('Data generation mode undefined');
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

    params.channels=val(find(var=='channels')); // Number of channels in template file
    params.recordingRate=val(find(var=='recordingRate')); // Sampling rate in [Hz]
    params.recordingLength=val(find(var=='recordingLength')); // Length of recording in [ms]
    params.spikeTypes=val(find(var=='spikeTypes')); // Number of spike types in template file
    params.spikeLength=val(find(var=='spikeLength')); // Length of a spike in [samples]
    params.spikesPerType=val(find(var=='spikesPerType')); // Number of spikes to generate for each type
    params.maxAmplitude=val(find(var=='maxAmplitude')); // If not NaN, spike primary amplitudes normalized to this value [uV]
    params.noiseAmplitude=val(find(var=='noiseAmplitude')); // Standard deviation of Gaussian noise [uV]
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


// Normalize cluster labels found during spike sorting
function newclu=norm_clusters(aclu)

  u=unique(aclu);
  newclu=zeros(size(aclu,1),1);

  for i=1:size(u,1)
    newclu(find(aclu==u(i,1)))=i;
  end

endfunction


// Constructs the confusion matrix for a single clustering
function c=confusion(params,gpos,gclu,apos,aclu)

  aclu=norm_clusters(aclu); // Normalize cluster labels

  M=params.spikeTypes; // Number of original clusters
  N=max(aclu); // Number of clusters found

  c=zeros(M+1,N+1); // Define confusion matrix

  // Generated positions mark the start of a spike, obtained positions are
  // the points of threshold crossing, need to correct for offset
  offset=floor(params.spikeLength/2)-1;
  apos=apos-offset;

  // Fills out c(i,j) for i>0,j>0
  for i=1:M
    gp=gpos(find(members(gclu,i)),1);
    for j=1:N
      ap=apos(find(members(aclu,j)));
      pos=repmat(ap,1,2*params.tolerance+1)+repmat(-params.tolerance:1:params.tolerance,size(ap,1),1);
      tpos=gsort(pos(:),'g','i');
      c(i+1,j+1)=sum(members(gp,tpos));
    end
  end

  // Fills out c(i,0), i>0
  // For each original cluster i, finds all timestamps in apos that are from i, 
  // then calculates the size of the set difference between i and the acquired
  // elements
  for i=1:M
    gp=gpos(find(members(gclu,i)),1);
    pos=repmat(apos,1,2*params.tolerance+1)+repmat(-params.tolerance:1:params.tolerance,size(apos,1),1);
    tpos=gsort(pos(:),'g','i');
    c(i+1,1)=size(setdiff(gp,gp(find(members(gp,tpos)))),1);
  end

  //for i=1:M
  //    c(i+1,1)=params.spikesPerType-sum(c(i+1,:),2);
  //end

  // Fills out c(0,j), j>0
  for j=1:N
    c(1,j+1)=size(apos(find(members(aclu,j))),1)-sum(c(:,j+1),1);
  end

endfunction


// Returns the Completeness of each obtained cluster
function com=completeness(conf)

  k=params.penalty+1;

  for j=2:size(conf,2)
    com(1,j-1)=max(conf(2:$,j)./sum(conf(k+1:$,:),2));
  end

endfunction


// Returns the Purity of each obtained cluster
function pur=purity(conf)

  k=params.penalty+1;

  for j=2:size(conf,2)
    pur(1,j-1)=max(conf(2:$,j)./sum(conf(k+1:$,j),1));
  end

endfunction


//The function generates m element subsets of set 1:n
//Source: https://groups.google.com/forum/#!msg/comp.soft-sys.math.scilab/HyYwDyGeT5U/A_0u4YjDNlEJ
//Available in Matlab as nchoosek(1:n,m)
//Available in Scilab 6.0.x only with atomsInstall("specfun"), specfun_subset(1:n,m)
function [A] = subsets(n,m)
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

