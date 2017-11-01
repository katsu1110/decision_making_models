function lookinto_history

% without history term
yout0 = Yates_simple;
close all;
disp('Yates_smple.m was run')

% with history term
yout1 = Yates_history;
close all;
disp('Yates_history.m was run')

% visualize comparison
subplot(2,2,1)
plot(yout0.psth.stm_null, '-r')
hold on;
plot(yout0.psth.stm_pref, '-b')
ylabel('psth (w.o. history)')
title('stimulus')

subplot(2,2,2)
plot(yout0.psth.ch_null, '-r')
hold on;
plot(yout0.psth.ch_pref, '-b')
title('choice')

subplot(2,2,3)
plot(yout1.psth.stm_null, '-r')
hold on;
plot(yout1.psth.stm_pref, '-b')
ylabel('psth (wi. history)')

subplot(2,2,4)
plot(yout1.psth.ch_null, '-r')
hold on;
plot(yout1.psth.ch_pref, '-b')