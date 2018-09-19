function out=RippleDetection(RecordingValues, SleepScoring, ScoringArtefacts, Label, fsample)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spindle Detection Algorith for EEG/LFP siganls
% INPUT:    RecordingValues: Matrix containing Raw EEG or LFP values
%           SleepScoring: Matrix containing Sleep Scoring for each data point
%           ScoringArtefact: Matrix coding the artefacts (1) of the signal
%           for each data point.
%           Label: CHAR with the name of the Channel used.
%           fsample: sampling rate.
% OUT:      Output contains detection info in the following way:
%           out.label           = CELL containg the name of the Channel(s)
%           out.fsample         = sampling rate
%           out.trial           = Event(s) raw signal
%           out.time            = time/length of the Trial(s)
%           out.trialinfo       = Detailed info about each one of the trials.
%               out.trialinfo(:,1) = To which time bin corresponds that particular event 
%               out.trialinfo(:,2) = startTime: begin of event 
%               out.trialinfo(:,3) = midTime: center of event
%               out.trialinfo(:,4) = endTime: spindle: negative zero crossing
%               out.trialinfo(:,5) = duration: duration from start to end in seconds (spindle: betweet the two threshild crossings)
%               out.trialinfo(:,6) = maxTime: time of maximum (spindle: largest negative trough during spindle) in datapoints of original dataset
%               out.trialinfo(:,7) = minTime: time of minimum (spindle: largest positive peak during spindle) in datapoints of original dataset
%               out.trialinfo(:,8) = minAmp: amplitude of maximum (spindle: largest negative trough during spindle) in µV
%               out.trialinfo(:,9) = maxAmp: amplitude of minimum (spindle: largest positive peak during spindle) in µV
%               out.trialinfo(:,10)= p2pAmp: peak-to-peak amplitude (spindle: largest peak to largest trough) 
%               out.trialinfo(:,11)= p2pTime: time in seconds from peak-to-peak (spindle: abs(min-max)) 
%               out.trialinfo(:,12)= Power
%               out.trialinfo(:,13)= Frequency
%           out.vector          = Vector containing the presence of the events
%           out.sleeptimeinsec  = Time in SECONDS used for the Detections (excl. time spent in Artefacts)
%
% Authors:  Carlos N. Oyanedel - Uni Tübingen, Germany
%           Niels Niethard - Uni Tübingen, Germany
%           Thanks to Dr. Hong-Viet Ngo, University of Birmingham, UK
% Contact:  jan.born@uni-tuebingen.de
%
% Date:     17.09.2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Basic Parameters
% Ripple Duration in Seconds
MinRippleDur    = 0.025;
MaxRippleDur    = 0.5;
NumberofCycles  = 3;

% Filters in Hz
LowPassFilter   = 250;
HighPassFilter  = 150;

% Sleep Stage of Interest
%1: WAKE 2: NREM; 3:REM; 4:PREM 
SleepStages     = 1:4;
SleepStageToUse = 2;

% STD above mean
StdValue        = 2.5; %The amount of STD above the mean for the detection

% Trial Window in seconds
TrialWindow  = 2.5; % Window before and after the Min Peak

%% Filtering the signal
[d,e]               = butter(3,2*LowPassFilter/fsample,'low'); %Use butterworth filter 3rd order. 15Hz cutoff
FilteredLFP         = filtfilt(d,e,RecordingValues); %Filter Signal
[d,e]               = butter(3,2*HighPassFilter/fsample,'high'); %Use butterworth filter 3rd order. 8Hz cutoff
FilteredLFP         = filtfilt(d,e,FilteredLFP);
FilteredLFPEnvelope = abs(hilbert(FilteredLFP));
RecLength           = length(RecordingValues);

%% Smoothing the FilteredLFPEnvelope
run SmoothFilteredLFPEnvelope.m

%% Mean and STD of the SmoothedFilteredLFPEnvelope of an specific Sleep Stage
ValidSleepScoring                       = SleepScoring;
ValidSleepScoring(ScoringArtefacts==1)  = 0;
SmoothedFilteredLFPEnvelopeMean         = mean(SmoothedFilteredLFPEnvelope(find(ValidSleepScoring==SleepStageToUse)));
SmoothedFilteredLFPEnvelopeSTD          = std(SmoothedFilteredLFPEnvelope(find(ValidSleepScoring==SleepStageToUse)));

%% RippleVector
RippleVector    = zeros(length(RecordingValues),1);

for i=1:length(RecordingValues); %Detecting Ripples crossing the power threshold
    if SmoothedFilteredLFPEnvelope(i) > (SmoothedFilteredLFPEnvelopeMean+StdValue*SmoothedFilteredLFPEnvelopeSTD)
        RippleVector(i) = 1;
    end
end

%% If want to delete event detected in different sleep stages
% Delete detected values in the other Sleep Stages
SleepStagetoRemove                    = SleepStages(SleepStages~=SleepStageToUse);
RippleVector(ValidSleepScoring==SleepStagetoRemove(1)) = 0; %
RippleVector(ValidSleepScoring==SleepStagetoRemove(2)) = 0; %
RippleVector(ValidSleepScoring==SleepStagetoRemove(3)) = 0; %
RippleVector(ValidSleepScoring==8)    = 0; %In Artefacts 8
RippleVector(ScoringArtefacts==1)     = 0; %In Artefacts
RippleVector(ValidSleepScoring==0)    = 0; %In case of mistmatch in extending the Sleep Scoring, i.e. DiffScoringRec>0

%% Detecting Beginning and End of the ripple events
% Position of values detected
RippleVectorLoc = find(RippleVector==1);

% Beg and End of each event.
RippleEnd       = [];
RippleBeg       = [];

for i=2:length(RippleVectorLoc)-1
    if RippleVectorLoc(i) - RippleVectorLoc(i-1) > 1
        RippleBeg = [RippleBeg,RippleVectorLoc(i)];
    end
    if RippleVectorLoc(i+1) - RippleVectorLoc(i) > 1
        RippleEnd = [RippleEnd,RippleVectorLoc(i)];
        
    end
end
if size(RippleBeg,2) > size(RippleEnd,2)
    RippleBeg = RippleBeg(1:size(RippleEnd,2));
end
RippleBeg       = [RippleVectorLoc(1),RippleBeg];
RippleEnd       = [RippleEnd, RippleVectorLoc(end)];
RippleDetected  = [RippleBeg;RippleEnd];

%% Duration Threshold
% Looking for Ripples with a Min and Max length
tmpRippleDuration  = diff(RippleDetected);

MinRippleDur    = MinRippleDur*fsample;
MaxRippleDur    = MaxRippleDur*fsample;

j=1;
for i = 1:size(RippleDetected,2)
    if tmpRippleDuration(i) <= MaxRippleDur && tmpRippleDuration(i) >= MinRippleDur %Checking which ones fulfill the duration criteria
        RippleTmp1(1,j) = RippleDetected(1,i);
        RippleTmp1(2,j) = RippleDetected(2,i);
        j               = j+1;
    end
end

%% Cycle Threshold - at least 3 pos and 3 negative cycles.
% Positive Peaks
PosPeakRippleNumber = [];
for i = 1:length(RippleTmp1);
    PosPeakRippleNumber(:,i) = length(findpeaks(FilteredLFP(RippleTmp1(1,i):RippleTmp1(2,i))));
end

% Negative Peaks
NegPeakRippleNumber = [];
for i = 1:length(RippleTmp1);
    NegPeakRippleNumber(:,i) = length(findpeaks(-FilteredLFP(RippleTmp1(1,i):RippleTmp1(2,i))));
end

%% Valid Ripples
ValidRipples    = RippleTmp1(:,find(PosPeakRippleNumber >= NumberofCycles & NegPeakRippleNumber >= NumberofCycles)); % At least 3 positive peaks.
% Duration in sec
RippleDuration  = diff(ValidRipples)/1000;
% Positive Peaks
PeakRippleNumber = zeros(1,size(ValidRipples,2));
for i = 1:length(ValidRipples);
    PeakRippleNumber(:,i) = length(findpeaks(FilteredLFP(ValidRipples(1,i):ValidRipples(2,i))));
end
% Frequency
tmpFreq = PeakRippleNumber./RippleDuration;

%% Valid Ripple Vector
ValidRippleVector = zeros(1,size(ValidSleepScoring,1));
for i=1:size(ValidRipples,2)
    ValidRippleVector(1,ValidRipples(1,i):ValidRipples(2,i)) = 1;
end

%% Ripple Power
tmppower = zeros(size(ValidRipples,2),1);
for iPow = 1:size(ValidRipples,2);
    tmppower(iPow,1) = trapz(SmoothedFilteredLFPEnvelope(ValidRipples(1,iPow):ValidRipples(2,iPow)));
end

%% Detecting Local Min and Max
% Local Minima
MinValidRippleVal = []; % Local Minima Value Matrx
MinValidRippleLoc = []; % Local Minima Position
for i = 1:length(ValidRipples);
    [MinValidRippleVal(:,i), MinValidRippleLoc(:,i)] = min(FilteredLFP(ValidRipples(1,i):ValidRipples(2,i)));
end

% Local Maxima
MaxValidRippleVal = []; % Local Maxima Value Matrx
MaxValidRippleLoc = []; % Local Maxima Position
for i = 1:length(ValidRipples);
    [MaxValidRippleVal(:,i), MaxValidRippleLoc(:,i)] = max(FilteredLFP(ValidRipples(1,i):ValidRipples(2,i)));
end

% Positions related to the whole recording for the Min and Max for each Valid Ripple Detected
for i=1:size(ValidRipples,2);
    ValidRippleMinLoc(i) = (ValidRipples(1,i)+MinValidRippleLoc(i));
    ValidRippleMaxLoc(i) = (ValidRipples(1,i)+MaxValidRippleLoc(i));
end

% Valid Max2Min info (Time Difference)
Max2MinValidLoc     = [MinValidRippleLoc; MaxValidRippleLoc];
Max2MinValidDiffLoc = diff(Max2MinValidLoc);

% Valid Max2Min info (Amplitud)
Max2MinValidRipple      = [MinValidRippleVal; MaxValidRippleVal];
Max2MinValidRipplesDiff = diff(Max2MinValidRipple);

% Detecting Max Values and position for the Envelope of each Valid Ripple
MaxRippleEnvelopeVal = []; % Local Maxima Value Matrx for the Envelope
MaxRippleEnvelopeLoc = []; % Local Maxima Position for the Envelope
for i = 1:size(ValidRipples,2);
    [MaxRippleEnvelopeVal(:,i), MaxRippleEnvelopeLoc(:,i)] = max(SmoothedFilteredLFPEnvelope(ValidRipples(1,i):ValidRipples(2,i)));
end

% Positions related to the whole recording for the Envelope Valid Ripple Max Detected
for i=1:size(ValidRipples,2);
    ValidRippleEnvelopeMaxLoc(i) = (ValidRipples(1,i)+MaxRippleEnvelopeLoc(i));
end

%% Ripple Grand Average
RippleGrandAverage  = zeros(size(ValidRipples,2),5000);
out.trial           = {};
out.time            = {};
WindowGrandAverage  = TrialWindow*fsample; % Window: 2.5 sec before and after the Min Peak
Cut                 = 0;  %In case the Trial Window for the last detected event is longer than the recording 
CutBeg              = 0;  %In case the Trial Window for the first detected starts before the recording file
for i = 1:size(MinValidRippleLoc,2)
    if ((ValidRipples(1,i)+MinValidRippleLoc(i))+WindowGrandAverage-1) > size(RecordingValues,1)
        Cut = 1;
    elseif ((ValidRipples(1,i)+MinValidRippleLoc(i))-WindowGrandAverage-1) < 0
        CutBeg = 1;
    else
        RippleGrandAverage(i,:) = RecordingValues((ValidRipples(1,i)+MinValidRippleLoc(i))-WindowGrandAverage:(ValidRipples(1,i)+MinValidRippleLoc(i))+WindowGrandAverage-1);
        out.trial{1,i}          = RippleGrandAverage(i,:)*1000;
        out.time{1,i}           = -TrialWindow:1/fsample:TrialWindow-1/fsample;
    end
end

%% Creating eventdata
out.label       = {Label}; % Creating eventdata file
out.fsample     = fsample;
tmptrialinfo    = zeros(size(ValidRipples,2),13);
for iTrial=1:size(ValidRipples,2);
    tmptrialinfo(iTrial,1)   = (ceil(ValidRipples(1,iTrial)/1800000)*30); %To which timebin corresponds
    tmptrialinfo(iTrial,2)   = ValidRipples(1,iTrial); %startTime: begin of event (ripple: positive threshold crossing)
    tmptrialinfo(iTrial,3)   = (ValidRipples(1,iTrial)+ValidRipples(2,iTrial))/2; %midTime: center of event (ripple: largest trough)
    tmptrialinfo(iTrial,4)   = ValidRipples(2,iTrial); %endTime: (ripple: negative zero crossing) 
    tmptrialinfo(iTrial,5)   = RippleDuration(iTrial); %duration: duration from start to end in seconds (ripple: betweet the two threshild crossings)
    tmptrialinfo(iTrial,6)   = ValidRippleMinLoc(iTrial); %maxTime: time of maximum (ripple: largest negative trough during ripple) in datapoints of original dataset
    tmptrialinfo(iTrial,7)   = ValidRippleMaxLoc(iTrial); %minTime: time of minimum (ripple: largest positive peak during ripple) in datapoints of original dataset
    tmptrialinfo(iTrial,8)   = MinValidRippleVal(iTrial)*1000; %minAmp: amplitude of maximum (ripple: largest negative trough during ripple) in ï¿½V
    tmptrialinfo(iTrial,9)   = MaxValidRippleVal(iTrial)*1000;%maxAmp: amplitude of minimum (ripple: largest positive peak during ripple) in ï¿½V
    tmptrialinfo(iTrial,10)  = Max2MinValidRipplesDiff(iTrial)*1000; %p2pAmp: peak-to-peak amplitude (ripple: largest peak to largest trough)
    tmptrialinfo(iTrial,11)  = abs(Max2MinValidDiffLoc(iTrial))/fsample; %p2pTime: time in seconds from peak-to-peak (ripple: abs(min-max)) 
end
tmptrialinfo(:,12)  = tmppower; %Power of each event detected
tmptrialinfo(:,13)  = tmpFreq'; %Freq of each event detected

if Cut == 1
    tmptrialinfo = tmptrialinfo(1:size(out.time,2),:);
end

if CutBeg == 1
    out.trial   = out.trial(2:end);
    out.time    = out.time(2:end);
end

out.trialinfo       = tmptrialinfo;
out.vector          = ValidRippleVector;
out.sleeptimeinsec  = ((size(find(SleepScoring == SleepStageToUse),1)) - (size(find(SleepScoring == SleepStageToUse & ScoringArtefacts == 1),1)))/fsample;

end