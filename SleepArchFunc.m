function out = SleepArchFunc(ScoringFile, RecordingLength, fsample)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:    ScoringFile: Nx2 matrix containing Scoring (:,1) and Artefact (:,2) Info
%           RecordingLength: Length of the EEG/LFP File
%           fsample: sampling rate.
% OUT:      Several Basic SleepArch parameters
%
% Authors:  Carlos N. Oyanedel - Uni Tübingen, Germany
%           Niels Niethard - Uni Tübingen, Germany
%           Thanks to Dr. Hong-Viet Ngo, University of Birmingham, UK
% Contact:  jan.born@uni-tuebingen.de
%
% Date:     17.09.2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Getting 
RawScoring  = ScoringFile(:,2);
RawArtefact = ScoringFile(:,3);

DiffScoringRec  =  RecordingLength - (length(RawScoring)*fsample*10);

% Extending RawScoring to the length of the recording
SleepScoring = zeros((RecordingLength),1);
parfor i = 1:RecordingLength-DiffScoringRec
    SleepScoring(i,1) = RawScoring(ceil((1/fsample)*i/10));
end

% Defining a specific time for analysis

REM     = find (SleepScoring==3); %find REM episodes
PREM    = find (SleepScoring==4); %find PREM episodes
NREM    = find (SleepScoring==2); %find NREM episodes
WAKE    = find (SleepScoring==1); %find Wake episodes

% REM
if isempty(REM)
    REMEpisodes = [];
else
    REMEndEpisode = [];
    REMBegEpisode = [];
    
    parfor i=2:length(REM)-1
        if REM(i) - REM(i-1) > 1
            REMBegEpisode = [REMBegEpisode,REM(i)];
        end
        if REM(i+1) - REM(i) > 1
            REMEndEpisode = [REMEndEpisode,REM(i)];
            
        end
    end
    
    REMBegEpisode   = [REM(1),REMBegEpisode];
    REMEndEpisode   = [REMEndEpisode, REM(end)];
    
    REMEpisodes     = [REMBegEpisode;REMEndEpisode];
end

NumberOfREMEpochs = size(REMEpisodes,2);

SleepArchitecture.REMEpisodesDuration           = diff(REMEpisodes);
SleepArchitecture.REMEpisodesDurationSeconds    = diff(REMEpisodes)/fsample; %Duration in seconds

TotalTimeMinREM   = sum(SleepArchitecture.REMEpisodesDurationSeconds)/60;

% PREM
if isempty(PREM)
    PREMEpisodes = [];
else
    
    PREMEndEpisode = [];
    PREMBegEpisode = [];
    
    parfor i=2:length (PREM)-1
        
        if PREM(i) - PREM(i-1) > 1
            PREMBegEpisode = [PREMBegEpisode,PREM(i)];
        end
        if PREM(i+1) - PREM(i) > 1
            PREMEndEpisode = [PREMEndEpisode,PREM(i)];
            
        end
    end
    
    PREMBegEpisode  = [PREM(1),PREMBegEpisode];
    PREMEndEpisode  = [PREMEndEpisode, PREM(end)];
    
    PREMEpisodes    = [PREMBegEpisode;PREMEndEpisode];
    
end

NumberOfPREMEpochs = size(PREMEpisodes,2);

SleepArchitecture.PREMEpisodesDuration          = diff(PREMEpisodes);
SleepArchitecture.PREMEpisodesDurationSeconds   = diff(PREMEpisodes)/fsample; %Duration in seconds

TotalTimeMinPREM = sum(SleepArchitecture.PREMEpisodesDurationSeconds)/60;

% NREM

if isempty(NREM)
    NREMEpisodes = [];
else
    
    NREMEndEpisode = [];
    NREMBegEpisode = [];
    
    parfor i=2:length (NREM)-1
        
        if NREM(i) -NREM(i-1) > 1
            NREMBegEpisode = [NREMBegEpisode,NREM(i)];
        end
        if NREM(i+1) - NREM(i) > 1
            NREMEndEpisode = [NREMEndEpisode,NREM(i)];
            
        end
    end
    
    NREMBegEpisode  = [NREM(1),NREMBegEpisode];
    NREMEndEpisode  = [NREMEndEpisode, NREM(end)];
    
    NREMEpisodes    = [NREMBegEpisode;NREMEndEpisode];
    
end

NumberOfNREMEpochs = size(NREMEpisodes,2);

SleepArchitecture.NREMEpisodesDuration          = diff(NREMEpisodes);
SleepArchitecture.NREMEpisodesDurationSeconds   = diff(NREMEpisodes)/fsample; %Duration in seconds

TotalTimeMinNREM    = sum(SleepArchitecture.NREMEpisodesDurationSeconds)/60;

% WAKE

WAKEEndEpisode = [];
WAKEBegEpisode = [];

parfor i=2:length (WAKE)-1
    
    if WAKE(i) - WAKE(i-1) > 1
        WAKEBegEpisode = [WAKEBegEpisode,WAKE(i)];
    end
    if WAKE(i+1) - WAKE(i) > 1
        WAKEEndEpisode = [WAKEEndEpisode,WAKE(i)];
        
    end
end

WAKEBegEpisode  = [WAKE(1),WAKEBegEpisode];
WAKEEndEpisode  = [WAKEEndEpisode, WAKE(end)];

WAKEEpisodes    = [WAKEBegEpisode;WAKEEndEpisode];

NumberOfWAKEEpochs = size(WAKEEpisodes,2);

SleepArchitecture.WAKEEpisodesDuration          = diff(WAKEEpisodes);
SleepArchitecture.WAKEEpisodesDurationSeconds   = diff(WAKEEpisodes)/fsample; %Duration in seconds

TotalTimeMinWAKE = sum(SleepArchitecture.WAKEEpisodesDurationSeconds)/60;

%% Summary and outputs
% Defining Output
SleepStage  = {'Wake';'NREM';'PREM';'REM'};
out.fsample = fsample;

% Number of Episodes
out.NumberOfEpisodes        = SleepStage;
out.NumberOfEpisodes{1,2}   = NumberOfWAKEEpochs;
out.NumberOfEpisodes{2,2}   = NumberOfNREMEpochs;
out.NumberOfEpisodes{3,2}   = NumberOfPREMEpochs;
out.NumberOfEpisodes{4,2}   = NumberOfREMEpochs;

% Onset
out.OnsetMin        = SleepStage;
out.OnsetMin{1,2}   = 0;
out.OnsetMin{2,2}   = (NREM(1,1)-1)/fsample/60;
if isempty(PREM)
else
    out.OnsetMin{3,2}   = (PREM(1,1)-1)/fsample/60;
end
if isempty(REM)
else
    out.OnsetMin{4,2}   = (REM(1,1)-1)/fsample/60;
end

% Mean Durantion +- SEM Duration
out.MeanDurationSec         = SleepStage;
out.MeanDurationSec{1,2}    = mean(SleepArchitecture.WAKEEpisodesDurationSeconds);
out.MeanDurationSec{2,2}    = mean(SleepArchitecture.NREMEpisodesDurationSeconds);
out.MeanDurationSec{3,2}    = mean(SleepArchitecture.PREMEpisodesDurationSeconds);
out.MeanDurationSec{4,2}    = mean(SleepArchitecture.REMEpisodesDurationSeconds);

% SEM for the Duration in Seconds
out.MeanDurationSec{1,3}    = std(SleepArchitecture.WAKEEpisodesDurationSeconds)/sqrt(length(SleepArchitecture.WAKEEpisodesDurationSeconds));
out.MeanDurationSec{2,3}    = std(SleepArchitecture.NREMEpisodesDurationSeconds)/sqrt(length(SleepArchitecture.NREMEpisodesDurationSeconds));
out.MeanDurationSec{3,3}    = std(SleepArchitecture.PREMEpisodesDurationSeconds)/sqrt(length(SleepArchitecture.PREMEpisodesDurationSeconds));
out.MeanDurationSec{4,3}    = std(SleepArchitecture.REMEpisodesDurationSeconds)/sqrt(length(SleepArchitecture.REMEpisodesDurationSeconds));

% Total Sleep Time in Minutes
out.TotalTimeMin        = SleepStage;
out.TotalTimeMin{1,2}   = TotalTimeMinWAKE;
out.TotalTimeMin{2,2}   = TotalTimeMinNREM;
out.TotalTimeMin{3,2}   = TotalTimeMinPREM;
out.TotalTimeMin{4,2}   = TotalTimeMinREM;

% Episodes
out.WAKEEpisodes    = WAKEEpisodes;
out.NREMEpisodes    = NREMEpisodes;
out.PREMEpisodes    = PREMEpisodes;
out.REMEpisodes     = REMEpisodes;

%% Scoring
out.SleepScoring                = [RawScoring,RawArtefact];