function Yates_neuron_stm(yout)

% preset params
len_frame = 1050;
nbin = 7;
frameperbin = len_frame/nbin;

% contrast
offset = 100;
co = 1*[zeros(1,offset),ones(1,len_frame),zeros(1,offset)];

% stimuli
s = [0 0.1 0.2 0.3 0.4 0.5];
for i = 1:length(s)
    para.adaptation(1).stm(i).seq = [zeros(1,offset), max(s)*ones(1,frameperbin*6), s(i)*ones(1, frameperbin), zeros(1,offset)];
    para.adaptation(2).stm(i).seq = [zeros(1,offset), min(s)*ones(1,frameperbin*6), s(i)*ones(1, frameperbin), zeros(1,offset)];
end
para.pref.stm  = [zeros(1, offset), 0.1*ones(1, frameperbin*7), zeros(1,offset)];
para.null.stm  = [zeros(1, offset), -0.1*ones(1, frameperbin*7), zeros(1,offset)];

% neural responses
rp = exp(conv(yout.kernel_stm, para.pref.stm) + conv(yout.kernel_co, co));
rn = exp(conv(yout.kernel_stm, para.null.stm) + conv(yout.kernel_co, co));

% visualize
time = [1:length(para.adaptation(1).stm(i).seq)] - offset;
close all;
figure;
subplot(1,3,1)
col = lines(length(s));
for i = 1:length(s)
    plot(time, para.adaptation(1).stm(i).seq, '-', 'color', 'k')
    hold on;
    plot(time, para.adaptation(2).stm(i).seq, '-', 'color', 'g')
    hold on;
end
title('stimulus')

subplot(1,3,2)
for i = 1:length(s)
    fr1 = exp(conv(yout.kernel_stm, para.adaptation(1).stm(i).seq)...
        + conv(yout.kernel_co, co));
    para.adaptation(1).stm(i).fr = fr1;
    fr2 = exp(conv(yout.kernel_stm, para.adaptation(2).stm(i).seq)...
        + conv(yout.kernel_co, co));
    para.adaptation(2).stm(i).fr = fr2;
    plot(time, fr1(1:length(time)), '-k')
    hold on;
    plot(time, fr2(1:length(time)), '-g')
    hold on;
end
title('psth')
legend('strong adap.', 'weak adap.')
legend('boxoff')

subplot(1,3,3)
for i = 1:length(s)
    plot(s(i), mean(para.adaptation(1).stm(i).fr(offset+len_frame-frameperbin:offset+len_frame)), 'ok')
    hold on;
    plot(s(i), mean(para.adaptation(2).stm(i).fr(offset+len_frame-frameperbin:offset+len_frame)), 'og')
    hold on;
end
title('tuning curve')

% subplot(2,2,4)
% plot(time, rp(1:length(time)), '-b')
% hold on;
% plot(time, rn(1:length(time)), '-r')
% legend('pref stm', 'null stm')
% legend('boxoff')