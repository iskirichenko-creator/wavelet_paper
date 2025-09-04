%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initial Conditions and Definitions

signal1=dt_i;                            % Analyzed signal
X=3.8415; % Chi-square quantile value 95
dt=1;                                    % Discretization step 
pseudo=1.23;                             % For Morlet wavelet, non-complex
                                         % pseudo=1.05; % For complex wavelet
xstart=0;
k=0;                                     % Number of samples from the start if not starting from the beginning
x0=xstart+k*dt;
V=1;             % Sedimentation rate 
alfa1 = autocorr(signal1,2);              % First two terms of the autocorrelation series
alfa=(alfa1(2,1)+alfa1(3,1)^0.5)/2;
Disp=var(signal1);                         % Variance of the series

n=length(signal1);
xlim = [0,n-1]*dt+x0;
Timelime=([0,n-1]*dt+x0)/V;
Time=((0:1:n-1)*dt+x0)/V;
time1=2012-((0:1:n-1)*dt+x0)/V;
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cone of Influence
a_max=n/(2*2^0.5);                         % Maximum wavelet decomposition scale
%a_max=150;                               % Manual setting of maximum scale
a_min=1;
delta=10;                                   % Wavelet decomposition step
b_1=((0:1:n/2)*dt+x0);                     % Abscissa of the left COI branch
b_2=(0:1:n/2)*dt;                          % Abscissa in the left COI branch
b1=(((n/2)+1:1:n)*dt+x0);                  % Abscissa of the right COI branch
b1_1=((n/2)+1:1:n)*dt;                     % Abscissa in the right COI branch
b3=(0:1:n-1)*dt+x0;                        % General abscissa
b3_1=((0:1:n-1)*dt+x0)/V;
COI=((b_2)/(2^0.5)).*(pseudo/V); COI1=((n*dt-b1_1)./(2^0.5)).*(pseudo/V);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cone of Influence 


scales=(a_min:delta:a_max);                   % Scale
coefs=cwt(dt_i,scales,'morl');       % Wavelet decomposition coefficients

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Normalization to the maximum and minimum term in each decomposition row
rowmin = min(coefs,[],2);
rowmax = max(coefs,[],2);
Brow = rescale(coefs,'InputMin',rowmin,'InputMax',rowmax);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Normalization to the maximum and minimum term in each decomposition row

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Confidence Levels
w=(abs(coefs)).^2;                            % Power spectrum of wavelet coefficients
w1=2*w./(X*Disp);                             % Power spectrum of wavelet coefficients
period = scales.*1.23;                        % Relationship between wavelet frequency and period
period1=period';
period2=period1.*dt/V;
freq = 1./ period;
p = (1-alfa^2) ./ (1-2*alfa*cos(freq*2*pi)+alfa^2); % Red noise model
p=p';
Spectr_teor = ((p)*(ones(1,n)));
Spectr_teor1 = w1./Spectr_teor;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Without normalization + Confidence Levels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fourier Analysis
T=1/n;
t=(0:n-1)*T;
NFFT = 2^nextpow2(n); % Next power of 2 from length of y
Y = fft(signal1,NFFT)/n;
Y=(abs(Y).^2)*(n/Disp);
Y1= fft(signal1,NFFT)/n;
Y1=(abs(Y1).^2)*(n/Disp);
f = linspace(0,0.5,NFFT/2+1);
f_p = dt./(f*V);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fourier Analysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plots

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Signal Plot
set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
ht = text(0.3,0.3,'Simple text','BackgroundColor','y');
set(ht,'Margin',20);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Signal Plot


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Global Wavelet Spectrum Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

global_wavelet=sum(w')/(n*Disp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DRAW SIGNAL
subplot('position',[0.1 0.70 0.65 0.25]);
plot(Time, i,'Color','k','LineWidth',1);
set(gca,'xLim',Timelime(:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlabel('Years');
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DRAW Wavelet Coefficient
subplot('position',[0.1 0.1 0.65 0.5]);
contour ( b3,period2,Brow,100, 'linewidth',0.5);
colormap(jet);
xlabel('Years');
ylabel('Period, Years');
set(gca, 'YScale', 'log') %log
hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DRAW Wavelet Coefficient

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DRAW Confidence Level
contour (b3,period2,Spectr_teor1,[-99,1],'k', 'linewidth',1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DRAW Confidence Level

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DRAW COI
plot(b_1,COI,'y','linewidth',2);
plot(b1,COI1,'y','linewidth',2);
ax = gca; % current axes
ax.YLim = [min(period2) max(period2)];
title('b) Wavelet Spectrum')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DRAW COI


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot Global Wavelet Spectrum
subplot('position',[0.77 0.1 0.2 0.5]);
plot(global_wavelet,period2, 'k', 'LineWidth',1);
ylim = [min(period2),max(period2)];
set(gca,'YLim',ylim(:))
hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot Global Wavelet Spectrum

plot((p*X),period2, 'Color','k','LineStyle','--');
ylim = [min(period2),max(period2)];
set(gca,'YLim',ylim(:))
hold on;
plot( Y1(1:NFFT/2+1),f_p,'Color','r','LineWidth',1);
plot(p,period2, 'Color','b','LineStyle','--');
set(gca, 'YScale', 'log')% log
hold on;
ylim = [min(period2),max(period2)];
set(gca,'YLim',ylim(:))
title('c) Fourier Power Spectrum')
hold off