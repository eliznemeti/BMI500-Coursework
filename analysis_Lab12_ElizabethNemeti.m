% reset the workspace
clear
close all

% load spiral drawing data
d = read_trc("lue-spiral.trc");

% set plotting parameters
TL = [0 5];
nr = 2;
nc = 3;

% plot the left hand marker in x-y-z
marker_name = "L.Finger3.M3";
marker_xyz = d{:,find(names(d) == "L.Finger3.M3") + [0:2]};

t = d{:,"Time"};
t_inds = t>min(TL)&t<max(TL);
t_secs = rem(t(t_inds),1)==0;

% first plot
figure
subplot(nr,nc,1)
hold on
plot(t(1:601), marker_xyz(1:601,1), ".-")
plot(t(1:601), marker_xyz(1:601,2), ".-")
plot(t(1:601), marker_xyz(1:601,3), ".-")
legend(["X", "Y", "Z"]);
xlim([0 5])
xlabel("seconds")
ylabel("mm")
title("Raw Data")
hold off

% second plot
subplot(nr,nc,2)
hold on
plot(marker_xyz(1:601,2), marker_xyz(1:601,3), "k")
xlabel("Y")
ylabel("Z")
title("Front View")
hold off

% Filter out large, slow movements with a high-pass butterworth filter at 2
% Hz cutoff and filter out jitter with a low-pass butterworth filter at 20
% Hz cutoff. A 6th order filter is fine.

% sampling freq fs is the reciprocal of the difference between two points
fs = 1/mean(diff(t));

% cutoff frequencies for the filter
fc_hi = 2/(fs/2);
fc_lo = 20/(fs/2);

% [b,a] = butter(n,Wn) returns the transfer function coefficients of an 
% nth-order lowpass digital Butterworth filter with normalized 
% cutoff frequency Wn [https://www.mathworks.com/help/signal/ref/butter.html]

%%% YOUR CODE HERE
[b_h, a_h] = butter(6, fc_hi); % creating butterworth filter
data_out_low = filtfilt(b_h, a_h, marker_xyz); % use filtfilt to filter data through

% third plot
subplot(nr,nc,3)
hold on
plot(data_out_low(1:601,2), data_out_low(1:601,3), "k") % to 601 so that it's up to 5 seconds
xlabel("Y")
ylabel("Z")
title("Low Frequency Component")
hold off

% calculate the first PC

%%% YOUR CODE HERE
[b_l, a_l] = butter(6, [fc_hi fc_lo]); % using butter filter at 6th order to make bandpass
data_out_band = filtfilt(b_l, a_l, marker_xyz); % using filtfilt on bandpass data
[coeff,score,latent] = pca(data_out_band); % calculating pca
Xcentered = score(:,1)*coeff(:,1)'; % centering our x points

% fourth plot
subplot(nr,nc,4)
hold on
plot(data_out_band(1:601,2), data_out_band(1:601, 3), '.k', 'MarkerSize',20)
p=polyfit(Xcentered(:,2),Xcentered(:,3),1); %  makes the line longer using y = mx +b 
x_pnt = linspace(-75,75, 1000); % x_pnt creates points to plot the line
y_pnt = polyval(p, x_pnt); % getting our y points from the equation
plot(x_pnt,y_pnt);
xlim([-75 75])
ylim([-75 75])
xlabel("Y")
ylabel("Z")
title("High Frequency Component and 1st PC")

% calculate projection onto first PC

%%% YOUR CODE HERE
proj = data_out_band*coeff(:,1); % projection onto 5th plot

% smooth with a savitsky-golay smoother
proj_smooth = smoothdata(proj,'sgolay');

% count zero crossings
zcd = dsp.ZeroCrossingDetector();
numZeroCross = cast(zcd(proj_smooth(t_inds)),"double");
tremorFrequency = (numZeroCross/2)/max(TL);

% get envelope from 25 sample moving average
env_width = 25;
env = movmax(proj_smooth(t_inds),env_width);

% use the median of the moving maximum as the estimator of the amplitude
amp = median(env);

ttl = round(tremorFrequency,1) + " Hz, " + round(2*amp,1) + " mm amplitude";

% plot
subplot(nr,nc,[5 6])
hold on
plot(t,proj,'k.')
plot(t,proj_smooth,'r')
h1 = refline(0,amp);
h2 = refline(0,-amp);
h1.Color = 0.5*[1 1 1];
h2.Color = 0.5*[1 1 1];
xlim(TL)
title(ttl)
ylabel("mm")
xlabel("seconds")
