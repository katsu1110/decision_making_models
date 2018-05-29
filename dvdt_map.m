function dvdt_map(res)
%% 
% simulate how decision variable and decision time in a race-model
% relates to confidence 
%
% INPUT: output struct from 'SDT_PKA.m'

close all;
h = figure;

% select a particular signal strength
idx = abs(res.assigned_stm)==res.stmMean;
assigned_stm = res.assigned_stm(idx);
category = res.category(idx);
accuracy = res.accuracy(idx);
decisiontime = res.decisiontime(idx);
dv1 = res.dv{1}(idx,:);
dv2 = res.dv{2}(idx,:);

% prepare for mapping between dv, dt, and probability correct
nstep_dv = 5; nstep_dt = 5;
% define integer grid of coordinates for the above data
[X,Y] = meshgrid(1:nstep_dv, 1:nstep_dt);
% define a finer grid of points
smoothstep = 0.05;
[X2,Y2] = meshgrid(1:smoothstep:nstep_dv, 1:smoothstep:nstep_dt);
% preallocation
mat = {zeros(nstep_dv, nstep_dt), zeros(nstep_dv, nstep_dt), ...
    zeros(nstep_dv, nstep_dt), zeros(nstep_dv, nstep_dt)}; 
m = 1;

% check the anti-correlation between the two integrators
meddt = median(decisiontime);
dv1temp = dv1(:, meddt);
dv2temp = dv2(:, meddt);
ntr = length(dv1temp);
dv1prc = prctile(abs(dv1temp), 0:(100/nstep_dv):100);
dv2prc = prctile(abs(dv2temp), 0:(100/nstep_dv):100);
for i = 1:nstep_dv
    dv1range = abs(dv1temp) >= dv1prc(i) & abs(dv1temp) < dv1prc(i+1);
    for k = 1:nstep_dv 
        dv2range = abs(dv2temp) >= dv2prc(k) & abs(dv2temp) < dv2prc(k+1);
        mat{m}(i,k) = log(sum(dv1range==1 & dv2range==1)...
            /ntr);
    end
end
% smooth
mat{m} = interp2(X, Y, mat{m}, X2, Y2, 'linear');
    
% extract dv, dt with respect to choice correctness
ntr = size(category, 1);
dt = nan(2, ntr);
dv = nan(2, ntr);
loc = [2,1];
for n = 1:ntr
    l = loc(accuracy(n)+1);
    dt(l,n) = decisiontime(n);
    dv(l,n) = min([abs(dv1(n,dt(l,n))), abs(dv2(n,dt(l,n)))]);
end
for l = 1:2
    dtprc = prctile(dt(l,:), 0:(100/nstep_dt):100);
    dvprc = prctile(dv(l,:), 0:(100/nstep_dv):100);
    for i = 1:nstep_dv
        dvrange = dv(l,:) >= dvprc(i) & dv(l,:) < dvprc(i+1);
        for k = 1:nstep_dt 
            dtrange = dt(l,:) >= dtprc(k) & dt(l,:) < dtprc(k+1);
            mat{m+l}(i,k) = sum(dtrange==1 & dvrange==1)...
                /ntr;
        end
    end
    % smooth
    mat{m+l} = interp2(X, Y, mat{m+l}, X2, Y2, 'linear');
end
m = m + l + 1;
mat{m} = log(mat{m-2} + mat{m-1});
mat{m-2} = log(mat{m-2});
mat{m-1} = log(mat{m-1});
mat{m+1} = mat{m-2} - mat{m-1};
crange = [min([mat{m-2}(:); mat{m-1}(:)]), max([mat{m-2}(:); mat{m-1}(:)])];
tlab = {'log(P(v|t,S))', 'log(P(correct,v,t|v=B,S))', 'log(P(error,v,t|v=B,S))', 'log(P(v,t|S))','log odds correct'};
subplot(2,3,1)
plot(mean(dv1,1),'-');
hold on;
plot(mean(dv2,1)','-');
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
for n = 1:m+1
    subplot(2,3,n+1)
    imagesc(mat{n})
    colorbar('southoutside');
    if n==1
        xlabel('accumulator 1')
        ylabel('accumulator 2')
    else
        if ismember(n,[2,3])
            caxis(crange)
        end
        xlabel('decision time')
        if ismember(n,[2, 4])
            ylabel({'acculumated evidence','of losing accumulator'})
        end
    end
    title(tlab{n})
    set(gca, 'YDir', 'normal')
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
end