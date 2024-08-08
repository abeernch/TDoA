freq = 2.5e9;
c = 3e8;
lambda = c/freq;
pri = 667e-6;
pw = 3e-6;
np = 72;
samp_rate = 2e6;
samp_per_pw = round(samp_rate*pw);
samp_per_pri = round(samp_rate*pri);

tx_sig = repmat([ones(1,samp_per_pw) zeros(1,samp_per_pri-samp_per_pw)],1,np);
rx_sig = circshift(tx_sig,500);


t = 0:1/samp_rate:pri*np-1/samp_rate;
K = 2^(nextpow2(np));
dop_lim = lambda*(1/pri)/4;
%%
load('mti filter time domain 7th order.mat')
peak_detector_param.fastTime_peak_neighbour = 3;
peak_detector_param.slowTime_peak_neighbour = 3;
peak_detector_param.slowTime_target_Max_diff_threshold = 15;
peak_detector_param.initial_Blind_RangeBin = 20;
radar.lfm_sll.en = false;
r_ax = linspace(0,pri*150*1e3,samp_per_pri);
d_ax = linspace(-dop_lim,dop_lim,K);
%%
f = figure;ax = axes(f);
f2 = figure;ax2 = axes(f2);
%%
for v = 1:1:240
fd = (2*v)/lambda;
s = exp(2j*pi*fd*t);
sig_pow = db2mag(-100);
sig_final = awgn(tx_sig+sig_pow*rx_sig.*s,110);
sig_2d = reshape(sig_final,samp_per_pri,np);
mf_2d = conv2(sig_2d,ones(samp_per_pw,1),'full');mf_2d = mf_2d(samp_per_pw:end,:);
mti_2d = filter(Num,1,mf_2d,[],2);
wind =transpose(chebwin(np,70)*ones(1,samp_per_pri));
fft_2d_mf = db(abs(fftshift(fft(mf_2d.*wind,K,2),2)));
fft_2d = fftshift(fft(mti_2d.*wind,K,2),2);
[rpeaks,ind_x] = max(db(abs(fft_2d)),[],2);
[hcfar,hnoise,~,~] = rcvr_cfar_sos(rpeaks,50,16,0.75,8);
[targets,~] =  rcvr_peak_detector_Generic_mianwali_2021(db(abs(fft_2d)),rpeaks.', hcfar.',ind_x.',peak_detector_param,fft_2d_mf.',radar,45);

% CFAR 

plot(ax,r_ax,[rpeaks hcfar],'LineWidth',1.2),ax.NextPlot = 'add';
if(~isempty(targets))
plot(ax,r_ax(targets(1,:)),rpeaks(targets(1,:)),'ro','LineWidth',2)
end
ax.NextPlot = 'replace';
grid on,grid minor
ylim(ax,[-90 -50]),xlim(ax,[0 r_ax(end)])
title(ax,'CFAR Target Detection','FontSize',20)
% FFT TARGET

mesh(ax2,d_ax,r_ax,db(abs(fft_2d)))
zlim(ax2,[-180 -50]),clim(ax2,[-120 -50])
axis(ax2,'tight')
view(ax2,[90 90])
% if(~isempty(targets))
% plot(ax,r_ax(targets(1,:)),rpeaks(targets(1,:)),'ro','LineWidth',2)
% end
% grid on,grid minor
% ylim([-90 -50]),xlim([0 r_ax(end)])
% title('CFAR Target Detection','FontSize',20)
pause(0.2)
end