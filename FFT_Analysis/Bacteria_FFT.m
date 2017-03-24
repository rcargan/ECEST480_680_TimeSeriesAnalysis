function Bacteria_FFT( bact_name, StoolA, SalivaA, range, dir, save_plots)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Fill in Missing Values

if verLessThan('matlab','9.1')
    % Replace NaN Values with 0s
    StoolA_NaNs = find(isnan(StoolA));
    StoolA(StoolA_NaNs) = 0;
    SalivaA_NaNs = find(isnan(SalivaA));
    SalivaA(SalivaA_NaNs) = 0;
else
    % Linearly Interpolate Values for Missing Stool Samples
    StoolA = fillmissing(StoolA,'linear');
    % Replace Missing Saliva Samples with 0s
    SalivaA = fillmissing(SalivaA,'constant',0);
end

%% Figure Formatting
range_text = [num2str(range(1)), '_', num2str(range(end))];
figure('name',bact_name);
%FigHandle = figure('name',bact_name,'Position', [100, 100, 300, 600]);
%hold on;


%% Plot Time Series
subplot(2,2,1);
plot(range,StoolA);
title(['Stool A Abundance, Day ', range_text],'FontSize',16);
xlabel('Time (days)','FontSize',14);
ylabel('Occurances','FontSize',14);

subplot(2,2,2);
plot(range,SalivaA,'g');
title(['Saliva A Abundance, Day ', range_text],'FontSize',16);
xlabel('Time (days)','FontSize',14);
ylabel('Occurances','FontSize',14);


%% Calculate and Plot FFTs
% Find Length of Total Signal
T = length(StoolA);
% Make Number of Elements Odd for Cleaner Processing
if(~mod(T,2))
    StoolA = StoolA(1:end-1);
    T = T-1;
end

% Create a Linearly Spaced Frequency Vector & Convert it to Period Length
fft_x = linspace(-T/2,T/2,T);
fft_x = T./fft_x;

% Calculate FFT to Account for Phase Shift, Remove DC Component, and
% Pull Out 3 Largest Periods
fftA = abs(fftshift(fft(StoolA)));
fftA(ceil(T/2)) = 0;
[sorted, I] = sort(fftA,'descend');
Per = fft_x(I([2,4,6]));

% Plot FFT
subplot(2,2,3);
stem(fft_x,fftA,'o');
title(['Stool A Periods = ',num2str(Per(1),3),', ',num2str(Per(2),3),', ',num2str(Per(3),3)],'FontSize',16);
xlabel('Period (days)','FontSize',14);
ylabel('Magnitude','FontSize',14);
axis([-65 65 0 inf])

% Find Length of Total Signal
T = length(SalivaA);

% Make Number of Elements Odd for Cleaner Processing
if(~mod(T,2))
    SalivaA = SalivaA(1:end-1);
    T = T-1;
end

% Create a Linearly Spaced Frequency Vector & Convert it to Period Length
fft_x = linspace(-T/2,T/2,T);
fft_x = T./fft_x;

% Calculate FFT to Account for Phase Shift, Remove DC Component, and
% Pull Out 3 Largest Periods
fftA_s = abs(fftshift(fft(SalivaA)));
fftA_s(ceil(T/2)) = 0;
[sorted, I] = sort(fftA_s,'descend');
Per = fft_x(I([2,4,6]));

% Plot FFT
subplot(2,2,4);
stem(fft_x,fftA_s,'go');
title(['Saliva A Periods = ',num2str(Per(1),3),', ',num2str(Per(2),3),', ',num2str(Per(3),3)],'FontSize',16);
xlabel('Period (days)','FontSize',14);
ylabel('Magnitude','FontSize',14);

%% Save Output Plots
if(save_plots)
    saveas(gcf,[dir, bact_name, range_text, '.png'])
end
    
end