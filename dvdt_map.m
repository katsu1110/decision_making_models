function dvdt_map(res)
%% 
% simulate how decision variable and decision time in a race-model
% relates to confidence 
%
% INPUT: output struct from 'SDT_PKA.m'

ntr = size(res.category, 1);
dt = nan(2, ntr);
dv = nan(2, ntr);
for n = 1:ntr
    if res.accuracy(n)
        l = 1;
    else
        l = 2;
    end
    dt(l,n) = res.decisiontime(n);
    if res.choice(n) == 1
        dv(l, n) = res.dv{2}(n,end);
    elseif res.choice(n) == 0
        dv(l, n) = res.dv{1}(n,end);
    end
end
close all;
h = figure;
subplot(1,2,1)
scatter(dt(1,:),abs(dv(1,:)),100,'s','markerfacecolor',[1 0 0],...
    'markerfacealpha',0.1,'markeredgecolor',[1 0 0],...
    'markeredgealpha',0.1)
xlabel('decision time')
ylabel({'acculumated evidence','of losing accumulator'})
title('correct trials')
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
subplot(1,2,2)
scatter(dt(2,:),abs(dv(2,:)),100,'s','markerfacecolor',[0 0 1],...
    'markerfacealpha',0.1,'markeredgecolor',[0 0 1],...
    'markeredgealpha',0.1)
xlabel('decision time')
ylabel({'acculumated evidence','of losing accumulator'})
title('error trials')
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')