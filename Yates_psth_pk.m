function [para] = Yates_psth_pk

% parameters for "grid search"
ntr = 200;
tmax = [150,300,400,500,600];
tau = [15,30,40,50,60];
kernel_stm = [0.01,0.02,0.04,0.08,0.1];
kernel_co = [0.01,0.02,0.04,0.08,0.1];

% loop for "grid search"
repeat = 50;
para = struct('repeat',[]);
label = [];
psth_stm = [];
psth_ch = [];
amp = [];
j = 1;
for i = 1:repeat
    disp(['iteration: ' num2str(i)])
    
    % randomly get index
    a = datasample(1:5,1);
    b = datasample(1:5,1);
    c = datasample(1:5,1);
    d = datasample(1:5,1);
    
    disp(['tmax: ',num2str(tmax(a)),', tau: ',num2str(tau(b)),...
        ', kernel_stm: ',num2str(kernel_stm(c)),', kernel_co: ',num2str(kernel_co(d))])
    try
        df = Yates_simple('ntr',ntr,'tmax',tmax(a),'tau',tau(b),'kernelgain_c',kernel_co(d),'kernelgain_s',kernel_stm(c),'nofig');
        
        label = [label, 1:7];
        psth_stm = [psth_stm, df.psth.stm_diff];
        psth_ch = [psth_ch, df.psth.ch_diff];
        amp = [amp, df.amp_diff];
        
        para.repeat(j).data = struct('df',df,'ntr',ntr,'tmax',tmax(a),'tau',tau(b),...
            'kernel_stm',kernel_stm(c),'kernel_co',kernel_co(d));
        disp('successfully run...')
        j = j + 1;
    catch
        disp('error occured...skipped.')
    end
    
end

% visualization
close all;
colors = lines(7);
figure;
for n = 1:7
    col = colors(n,:);
    subplot(1,2,1)
    scatter(psth_stm(label==n), amp(label==n), 10, 'filled', 'o', 'markerfacecolor', col, 'markeredgecolor',col,...
        'markerfacealpha',0.4,'markerfacealpha',0.8)
    xlabel('PSTH (stm: pref - anti)')
    ylabel('PK (conf: high - low)')
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
    hold on;
    
    subplot(1,2,2)
    scatter(psth_ch(label==n), amp(label==n), 10, 'filled', 'o', 'markerfacecolor', col, 'markeredgecolor',col,...
        'markerfacealpha',0.4,'markerfacealpha',0.8)
    xlabel('PSTH (ch: pref - anti)')
    ylabel('PK (conf: high - low)')
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
    hold on;
end
subplot(1,2,1)
xx = get(gca, 'XLim');
yy = get(gca, 'YLim');
plot([0 0],yy,':k')
hold on;
plot(xx,[0 0],':k')

subplot(1,2,2)
xx = get(gca, 'XLim');
yy = get(gca, 'YLim');
plot([0 0],yy,':k')
hold on;
plot(xx,[0 0],':k')

legend('1','2','3','4','5','6','7','location','eastoutside')
legend('boxoff')