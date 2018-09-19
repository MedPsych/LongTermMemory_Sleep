% Smoothing FilteredLFPEnvelope signal.
% MovingAverages = input(['Point window for smoothing the signal ']);%Point window for smoothing the signal
MovingAverages = 200;%Point window for smoothing the signal in ms
SmoothedFilteredLFPEnvelope = zeros(size(FilteredLFPEnvelope));

for i = 1:RecLength
    if i < MovingAverages/2 % SmoothedMean for the first 100 data points
        X = FilteredLFPEnvelope(1:(round(MovingAverages/2)+(i-1)));

        MeanData = mean (X);
        SmoothedFilteredLFPEnvelope (i) = MeanData;
    end
    
    if i >= MovingAverages/2 && i < RecLength - MovingAverages/2
        X = FilteredLFPEnvelope((i+1)-round(MovingAverages/2):round(MovingAverages/2)+(i-1));
 
        MeanData = mean (X);
        SmoothedFilteredLFPEnvelope (i) = MeanData;
    end
    if i >= RecLength - MovingAverages/2 % SmoothedMean for the last 100 data points
        X = FilteredLFPEnvelope((i+1)-round(MovingAverages/2):RecLength);
        
        MeanData = mean (X);
        SmoothedFilteredLFPEnvelope (i) = MeanData;
    end

end

clear MovingAverages FilteredLFPEnvelope i X MeanData