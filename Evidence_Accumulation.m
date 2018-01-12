function [dv, evidence] = Evidence_Accumulation(acc,e,ntime,ntrial,stim,model)

% function dv = Evidence_Accumulation(acc,ntime,ntrial)
% returns decision-variable for simulations of zero-signal case
% dv[number of time steps, number of acceleration parameters, number of
% trials]
% Allows one to obtain reaction time and choice for a particular bound B as
%      reaction time for first trial for first acc parameter: 
%               rt = find(dv(:,idx_acc,idx_trial)>2,1)
%               if isempty(rt), rt=size(dv,1); end
%      corresponding choice: 
%               choice = sign(dv(rt,idx_acc,idx_trial));

if nargin<1, acc=[0 0.01 -0.01]; end % [perfect, accelerate, deccelerate]
if nargin<2, e=1; end % average unsigned evidence per frame
if nargin<3, ntime=100; end % 100 time steps
if nargin<4, ntrial=1000; end % number of simulated trials
if nargin<5, stim='normal'; end % stimulus to be simulated (see below)
if nargin<6, model='linear'; end % model to be simulated (see below)

dv=zeros(ntime, length(acc),ntrial); % decision variable
switch stim
  case 'normal'
    evidence=randn(ntime,ntrial);
  case 'binary'
    evidence=-0.5+[ones(ntime/2,ntrial); zeros(ntime/2,ntrial)];
    for i=1:ntrial
      evidence(:,i)=evidence(randperm(ntime),i);
    end
  otherwise
    error(['invaid stim: ' stim]);
end
evidence=e*evidence;
evidence_upscale=ones(length(acc),1);
acc=acc'*ones(1,ntrial);
dv(1,:,:)=evidence_upscale*evidence(1,:); % starting point
for i=2:ntime
  x=reshape(dv(i-1,:,:),size(acc,1),ntrial);
  switch model
    case 'linear'
      x=x+acc.*x;       % dv=dv+acc*dv
    case 'sigmoid'
      x=x+acc.*tanh(x); % dv=dv+acc*tanh(dv)
    otherwise
      error(['invalid model: ' model]);
  end
  dv(i,:,:)=x+evidence_upscale*evidence(i,:);
end
