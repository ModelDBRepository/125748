function [spike_shapes, spike] = extract_spike(time, data, ddata, threshold, width, Tpre, Tpost, FLAG)
%   [SPIKE_SHAPES, SPIKE] = EXTRACT_SPIKE(TIME, DATA, DDATA, THRESHOLD, WIDTH, TPRE, TPOST, FLAG)
%
%     Version 2.0 - threshold-crossing triggered extraction of the peaks shape; no 'delaying' included.
%
%     Extract details on the occurrence time, the amplitude, the upward and downward speed and the peak shapes of a peak-train
%     by iterating a simple peak detection algorithm, based on (positive derivative) threshold-crossing.
%
%  *  The input argument  TIME refers to the data time column [s] and it is used for many internal operations and to extract the peak-times.
%  *  The input arguments DATA is the raw data trace [mV], while DDATA is the raw data trace to be used for detecting threshold-crossings.
%     In most of the situations DDATA is the time derivative of DATA (i.e. [mV/s], DDATA = diff(DATA,1)/dt, where dt = TIME(2)-TIME(1)). However it
%     can be DDATA = DATA or DDATA = zDATA, with zDATA a different raw data trace sharing the same length of DATA and TIME.
%  *  The input argument THRESHOLD, specified in the units of DDATA, is the (positive derivative) peak-event-detection threshold (e.g. 1.5E5 mV/s).
%  *  The input argument WIDTH is the minimal distance between two successive peak-events and it is assumed to be measured in [ms] (e.g. 2 ms).
%  *  The input arguments Tpre and Tpost specify, in [ms], the observation window to retrieve the shape of the detected spike.
%  *  The input argument FLAG specifies whether the detected spike shape to be put on SPIKE_SHAPES must be aligned according to the
%      threshold-crossing time (FLAG == 1), or alternatively to the time of the spike peak (FLAG == 0).
%
%  *  The output argument SPIKE_SHAPES is a cell array, containing the shape of each individual spike over a window specified by Tpre and Tpost.
%
%  *  The output argument SPIKE is a (N x 4) numeric array, containing respectively:
%                       - the absolute time of occurrence of the detected peak [s] (extracted from the corresponding element in TIME)
%                       - the amplitude of the peak [mV]                           (extracted from DATA)
%                       - the maximal upstroke slope of the peak [mV/s].           (extracted from the time derivative of DATA)
%                       - the maximal downstroke slope of the peak [mV/s].         (extracted from the time derivative of DATA)
%
%     � 2002 - Michele Giugliano, PhD (http://www.giugliano.info) (Bern, Friday July 5th, 2002 - 18:02)  (bug-reports to michele@giugliano.info)
%
if (isempty(time) | isempty(data) | isempty(ddata))
 disp(sprintf('Extract_Spike: (Error) empty imput data structures'));        
 return;
end % if
if (width <= 0)
 disp(sprintf('Extract_Spike: (Error) Wrong peak-event width size.'));        
 return;
end % if

spike_shapes = {};                              % The output data structure 'spike_shapes' is initialized.
spike     = [];                                 % The output data structure 'spike' is initialized. 
dt        = time(2) - time(1);                  % The sampling interval is evaluated (hp: it is expressed in [s]).
min_index = 1;        max_index = length(time); % Minimal and maximal sample index, respectively.

Iw        = round((width/1000.)/dt);            % Minimal number of samples between two successive peak-event (related to WIDTH).
Ipre      = round((abs(Tpre) /1000.)/dt);       % Number of samples, before the threshold crossing, to be extracted into 'spike_shapes{i}'.
Ipost     = round((abs(Tpost)/1000.)/dt);       % Number of samples, after  the threshold crossing, to be extracted into 'spike_shapes{i}'.

% Below, I took advantage of the MATLAB command 'find' to extract those part of the DDATA trace, whose amplitude is larger
% than an appropriately chosen threhsold, defining the DATA peak-events in terms of threshold-crossings of DDATA.

upstrokes                         = ddata;             % The vector 'upstrokes' is filled with the content of the DDATA signal.
upstrokes(find(ddata<=threshold)) = nan;               % Those elements in DDATA, whose amplitude is below the threshold, are marked as 'nan'.
Nup                               = length(upstrokes); % Here the length of 'upstrokes' is determined.

if (Nup == 0)
  disp(sprintf('Extract_Spike: (Error) No (positive derivative) threshold crossing detected!'));    
 return;                                               % If no (positive derivative) threshold crossing occurs, then the routine ends and returns an empty output.
end % if

if (upstrokes(1)~=nan)                                 % If the very first element of upstrokes is not a NaN
    upstrokes(1) = nan;                                % I conventionally set it to NaN so that the peak-event extraction algorithm will work anyway.
end % if

temp       = [];                                       % Temporary structure containing the temporary peak-event waveform..
tmp        = [];                                       % Temporary structure containing the temporary peak-event waveform..
counts     = 0;                                        % Counter of the peak-event number. 
t_last     = -inf;                                     % Last time a peak-event occurred [s] (i.e. it is initialized to a very remote time).
discard       = 0;                                     % Variable counting the number of times a peak has been discarded (just to trigger a warning comment only).

%---------------------------------------------------------------------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------------------------------------------------------------------
for J=2:Nup-1,                                         % I go through each element in 'upstrokes' to detect the (positive) threshold crossings. 
 if ( isnan(upstrokes(J-1)) & ~isnan(upstrokes(J)) )   % Here is the definition of "threshold crossing", occurring at time t*.
     if ((J-Iw)>=min_index) & ((J+Iw)<=max_index)      % If there are enough samples to fully include at least one peak-event then proceed.
     counts          = counts + 1;                     % Increase the current peak-event number counter.
     temp            = data(J-Iw:J+Iw);                % Extract from the raw data, samples between (t*-width) and (t*+width).
     hhh             = find(temp == max(temp));        % Extract the time (t*) at which the maximal value of the depolarization is reached ('the' peak-event).
     spike(counts,1) = time(J-Iw+hhh(1)-1);            % Write in the output data structure, the absolute time t* (note: index* = J-Iw+hhh(1)-1 ).
     spike(counts,2) = max(temp);                      % ..the amplitude of the peak-event (e.g. the depolarization at time t*),
     if ((J+hhh(1)-1-2*Iw)>=min_index) & ((J+hhh(1)-1)<=max_index)
      temp         = data(J+hhh(1)-1-2*Iw:J+hhh(1)-1);  % Extract from the raw data, (t* redefined to the peak time) the samples between t*-width and t*+width.s
      spike(counts,3) = max(diff(temp,1)/dt);           % the maximal positive derivative (max upstroke slope), within the specified window [t*-width ; t*+width]
      spike(counts,4) = min(diff(temp,1)/dt);           % the minimal derivative (max downstroke slope), within the specified window [t*-width ; t*+width]
      clear temp;                                       % Free some memory, as I don't need anymore the 'temp' data structure.
     
      if (FLAG)
       K = J;                                           % The index corresponding to the threshold-crossing time, to be used in the following..
      else
       K = J - Iw + hhh(1) - 1;                         % The index corresponding to the maximal amplitude of the spike is stored, to be used in the following..
      end %                                             % Once a spike has been detected, its shape must be extracted.

      if ((K-Ipre)<min_index) & ((K+Ipost)<=max_index)      % First, I must check whether I am running out the allowed range.                                                       
       spike_shapes{counts,1} = ((min_index-K):Ipost)*dt*1000.;                  
       spike_shapes{counts,2} = data(min_index:K+Ipost);    % I simply extract the raw waveform.
      elseif ((K-Ipre)<min_index) & ((K+Ipost)>max_index)        
       spike_shapes{counts,1} = ((min_index-K):max_index-K)*dt*1000.;                  
       spike_shapes{counts,2} = data(min_index:max_index);  % I simply extract the raw waveform.
      elseif ((K-Ipre)>=min_index) & ((K+Ipost)>max_index)        
       spike_shapes{counts,1} = (-Ipre:max_index-K)*dt*1000.;                  
       spike_shapes{counts,2} = data(K-Ipre:max_index);     % I simply extract the raw waveform.
      elseif ((K-Ipre)>=min_index) & ((K+Ipost)<=max_index)        
       spike_shapes{counts,1} = (-Ipre:Ipost)*dt*1000.;                  
       spike_shapes{counts,2} = data(K-Ipre:K+Ipost);       % I simply extract the raw waveform.
      end % if
     
      if ((spike(counts,1) - t_last) < (width/1000.))   % When the detected peak-event is unrealistically close to the preceding one, let's discard it.
       counts = counts - 1;                             % So let's decrement the peak-events counter.
       discard   = discard + 1;                         % A warning will be generated at the end of the routine.
      else                                              % Otherwise, when this is not the case,
       t_last = spike(counts,1);                        % simply keep trace of the present peak-event time for the next control, next cycle.
      end % if
     else                                               % (J+hhh(1)-1-2*Iw) < min_index or (J+hhh(1)-1) > max_index, so the algorithm cannot proceed in such a case.
      %disp(sprintf('Extract_Spike: (Warning) Not enough samples to proceed.'));
     end % if ((J+hhh(1)-1-2*I...   
    
    else                                               % (J-Iw) < min_index or (J+Iw) > max_index, so the algorithm cannot proceed in such a case.
     %disp(sprintf('Extract_Spike: (Warning) Not enough samples to proceed.'));
    end % if ((J-Iw)>...   
 end % if ( isnan(upstr..
end % for J=2:Nup-1
%---------------------------------------------------------------------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------------------------------------------------------------------
if (discard > 0) 
    disp(sprintf('Extract_Spike: (Warning) %d peaks have been discarded, being closer than %f.\n',discard,width)); 
end % if
disp(sprintf('Extract_Spike: %d peaks extracted.\n\n',counts));
