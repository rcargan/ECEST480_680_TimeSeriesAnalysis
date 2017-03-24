
%% Load in Data
Asa=readtable('DonorASaliva.csv'); %check to make sure 1001 columns present
allseqs=cell2mat((table2cell(Asa(1:end,4:1050))));  %convert all data to double
seq1=allseqs(:,1); %just look at 1st sequence
seq1=seq1(1:end-1);

%% Pairwise Distance

data=allseqs'; %transpose data so each row has gene sequence abundances
data=data(:,1:259);
data(isnan(data))=0; 

% normalize abundances
for i=1:size(data,1); %for each row
    %avg(i)= mean(data(i,:)); %compute average
    %stdev(i)= std(data(i,:)); %compute std
    maxval= max(data(i,:));
    normalized(i,:)= (data(i,:)/maxval);
    %normalized(i,:)= (data(i,:)-avg(i))/(stdev(i)); %store normalized data in new matrix
end

%signal filtering
peaks= zeros(size(normalized,1),1);
for i=1:size(normalized,1);
    peaks(i)= sum(normalized(i,:)>0.1); 
    pidx= find(peaks>10);     
end

distmat= pdist(normalized(1:400,4:end)); %take distances of each sequence
reshaped=squareform(distmat); %square matrix
figure
imagesc(reshaped); title('Gene Sequence Distance Matrix Visualization') %visualize

reshaped(reshaped==0)=inf; %turn all zeros to Inf

%find min value
for i=1:size(reshaped,1);
    [M(i),I(i)]=min(reshaped(i,:));
end

minvals= [M' I']; 
%small=find(minvals(:,1)<0.3);

    
%Visualize Differences To Choose Threshold 
subplot(3,1,1)
plot(normalized(154,:));
hold on
plot(normalized(191,:)+0.05);
title('Difference of 1.912','FontSize',14)
ylim([0 1]);
xlim([0 315]);
subplot(3,1,2);
plot(normalized(155,:));
hold on
plot(normalized(215,:)+0.05);
title('Difference of 0.824','FontSize',14)
ylim([0 1]);
xlim([0 315]);
subplot(3,1,3);
plot(normalized(169,:));
hold on
plot(normalized(149,:)+0.05);
title('Difference of 0.281','FontSize',14)
ylim([0 1]);
xlim([0 315]);

threshold= 0.7; %pick threshold

%create matrix of logicals
logicals= zeros(size(reshaped,1),size(reshaped,2));
for i=1:size(reshaped,1); %for each row
    for j=1:size(reshaped,2); %for each column
        if reshaped(i,j) < threshold
            logicals(i,j)=1;
        else 
            logicals(i,j)=0;
        end
    end
end

%find sequence with similarity to most other sequences
sim=zeros(size(logicals,1),1);
for i=1:size(logicals,1)
    sim(i)= sum(logicals(i,:));
end

row= logicals(245,:);
idx= find(row==1);

[ms,in]=max(sim);
simseqs=find(sim==ms) %sequences that have max similarity

for i=1:numel(simseqs); %for each sequence with high similarity with others...
    [~,col]=find(logicals(simseqs(i),:)==1); %find which sequences it is similar to
    figure
    for j=1:numel(col);
        plot(normalized(col(1,j)+j,:));
        %str=sprintf('Sequence #%d',simseqs(i));
        %title(str)
        hold on
        pause
     end
end







%% Conduct Autocorrelation and Differences
%seq1= data(294,:);
%seq1= data(343,:);
seq1=data(9,:);
[cr,lags]=xcorr(seq1,'coeff');
figure;

D1 = LagOp({1,-1},'Lags',[0,1]);
dY = filter(D1,seq1);
N=length(seq1);


subplot(2,2,1);
plot(seq1); xlim([0 N]); title('Sequence 1 Raw Data');
subplot(2,2,2);
stem(lags,cr,'.'); xlim([0 N]); title('Autocorrelation'); 
subplot(2,2,3);
plot(2:N,dY); xlim([0,N]); title('Data After One Difference');
subplot(2,2,4);
[cr2,lags2]=xcorr(dY,'coeff');
stem(lags2,cr2,'.');
title('Autocorrelation With One Difference'); 
xlim([0 N]);

% Mdl= arima('Constant',1,'D',0,'Seasonality',20,...
%  	'MALags',2,'SMALags',17, 'ARLags',4);
Mdl= arima('Constant',0,'D',0,'Seasonality',20,...
'MALags',3,'SMALags',20, 'ARLags',2, 'SARLags',1);
EstMdl = estimate(Mdl,seq1(130:240)');
%[Y,E] = simulate(EstMdl,314);
[yf, yMSE]= forecast(EstMdl,65,'Y0',seq1(130:200)');
figure
plot(seq1, 'Color',[0.75, 0.75, 0.75],'LineWidth',2);
hold on
h1=plot(2:1+65,yf,'r','LineWidth',2);
%h1= plot(250:249+65,yf, 'r','LineWidth',2); 
hold on
plot([70 70], [0 10000], '--k', 'LineWidth', 2)
plot([120 120], [0 10000], '--k', 'LineWidth', 2)
xlim([0 255])
ylim([0 10000])
legend('Original','Predicted')
title('ARIMA Prediction of Clostridiales Cluster','FontSize',13);
xlabel('Days','FontSize',13);
ylabel('Relative Abundance','FontSize',13)
dat= seq1(70:69+51);
