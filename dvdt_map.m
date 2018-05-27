function dvdt_map(res)
%% 
% simulate how decision variable and decision time in a race-model
% relates to confidence 
%
% INPUT: output struct from 'SDT_PKA.m'

% loop to extract dv, dt with respect to choice correctness
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

% map between dv, dt, and probability correct
nstep_dv = 5; nstep_dt = 5;

% define integer grid of coordinates for the above data
[X,Y] = meshgrid(1:nstep_dv, 1:nstep_dt);

% define a finer grid of points
smoothstep = 0.01;
[X2,Y2] = meshgrid(1:smoothstep:nstep_dv, 1:smoothstep:nstep_dt);

% map percent correct
mat = {zeros(nstep_dv, nstep_dt), zeros(nstep_dv, nstep_dt), zeros(nstep_dv, nstep_dt)}; 
for l = 1:2
    dtprc = prctile(dt(l,:), 0:(100/nstep_dt):100);
    dvprc = prctile(abs(dv(l,:)), 0:(100/nstep_dv):100);
    for i = 1:nstep_dv
        dvrange = abs(dv(l,:)) >= dvprc(i) & abs(dv(l,:)) < dvprc(i+1);
        for k = 1:nstep_dt 
            dtrange = dt(l,:) >= dtprc(k) & dt(l,:) < dtprc(k+1);
            mat{l}(i,k) = sum(dtrange==1 & dvrange==1)/ntr;
        end
    end
    % smooth
    mat{l} = interp2(X, Y, mat{l}, X2, Y2, 'linear');
end
mat{3} = mat{1}./(mat{1} + mat{2});

close all;
h = figure;
tlab = {'correct', 'error', 'odds'};
for n = 1:3
    subplot(1,3,n)
    imagesc(mat{n})
    c = colorbar('southoutside');
    if n==3
        c.Label.String = 'probability correct';
    else
        c.Label.String = 'probability';
    end
    xlabel('decision time')
    if n == 1
        ylabel({'acculumated evidence','of losing accumulator'})
    end
    title(tlab{n})
    xx = get(gca, 'XLim');
    yy = get(gca, 'YLim');
    set(gca, 'XTick', xx, 'XTickLabel', [0, round(max(dt(:)))])
    set(gca, 'YTick', yy, 'YTickLabel', [0, round(max(abs(dv(:))))])
    set(gca, 'YDir', 'normal')
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
end