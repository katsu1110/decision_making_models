function Yates_psth_pk

% parameters for "grid search"
ntr = 1000;
tmax = 600;
tau = 40;
kernel_stm =  [0.012,0.013,0.014,0.015,0.016,0.017,0.018];
kernel_co = [0.012,0.013,0.014,0.015,0.016,0.017,0.018];
offset_gain = [0.6, 0.7];
offset = 100;
len_frame = 1050;
nbin = 7;
frameperbin = len_frame/nbin;
pref  = [zeros(1, offset), 0.2*ones(1, frameperbin*7), zeros(1,offset)];
null  = [zeros(1, offset), -0.2*ones(1, frameperbin*7), zeros(1,offset)];
co = 1*[zeros(1,offset),ones(1,len_frame),zeros(1,offset)];

% yellow and green
y = [0.9576    0.7285    0.2285];
g = [0.1059    0.4706    0.2157];

% loop for "grid search"
time = [1:1250] - offset;
repeat = length(kernel_co);
close all;
for k = 1:2
    figure;
    for i = 1:repeat           
%         try
            if k==1
                coker = kernel_co(i);
                stker = 0.015;
            else
                coker = 0.015;
                stker = kernel_stm(i);
            end

            yout = Yates_simple_example('ntr',ntr,'tmax',tmax,'tau',tau,...
                'kernelgain_c',coker,'kernelgain_s',stker,...
                'offset_gain',offset_gain(1),'stm_gain',offset_gain(1),'nofig');
            disp(['tmax: ',num2str(tmax),', tau: ',num2str(tau),...
            ', kernel_stm: ',num2str(stker),', kernel_co: ',num2str(coker)])
        
            % visualize ==================
            % kernels
            subplot(7,repeat, i)
            plot(1:length(yout.kernel_stm), yout.kernel_stm, '-b')
            hold on;
            plot(1:length(yout.kernel_stm), -yout.kernel_stm, '-r')
            hold on;
            plot(1:length(yout.kernel_co), yout.kernel_co, '-k')
            xlim([1 length(yout.kernel_co)])
            set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
            if k==1
                title(coker)
            else
                title(stker)
            end
            if i==1
                ylabel('kernels')
            else
                set(gca, 'XTick', [])
            end
            
            % psth for stimulus (before nonlinearity)
            subplot(7,repeat,i+repeat)
            cv1 = conv(yout.kernel_co, co) + conv(yout.kernel_stm, pref);
            cv2 = conv(yout.kernel_co, co) + conv(yout.kernel_stm, null);
            plot(time, cv2(1:length(time)), '-r')
            hold on;
            plot(time, cv1(1:length(time)), '-b')
            xlim([1 length(time)])
            yy = get(gca, 'YLim');
            ylim([0 yy(2)])
            set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
            if i==1
                ylabel('ker * stm')
            else
                set(gca, 'XTick', [])
            end
            
            % psth for stimulus (after nonlinearity)
            subplot(7,repeat,i+2*repeat)
            fr1 = exp(cv1(1:length(time)));
            fr2 = exp(cv2(1:length(time)));
            plot(time, fr2, '-r')
            hold on;
            plot(time, fr1, '-b')
            xlim([1 length(time)])
            set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
            if i==1
                ylabel('FR')
            else
                set(gca, 'XTick', [])
            end
            
            % psth for stimulus (difference)
            subplot(7,repeat,i+3*repeat)
            plot(time, fr1 - fr2, '-k')
            xlim([1 length(time)])
            yy = get(gca, 'YLim');
            ylim([0 yy(2)])
            set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
            if i==1
                mes = sprintf('FR \n (pref - null)');
                ylabel(mes)
            else
                set(gca, 'XTick', [])
            end
            
            % psth for choice
            subplot(7,repeat,i+4*repeat)
            plot(time, yout.psth.ch_null(1:length(time)), '-r')
            hold on;
            plot(time, yout.psth.ch_pref(1:length(time)), '-b')
            xlim([1 length(time)])
            set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
            if i==1
                mes = sprintf('psth \n (choice)');
                ylabel(mes)
            else
                set(gca, 'XTick', [])
            end
            
            % psth for choice
            subplot(7,repeat,i+5*repeat)
            plot(time, yout.psth.ch_pref(1:length(time)) - ...
                yout.psth.ch_null(1:length(time)), '-k')
            xlim([1 length(time)])
            yy = get(gca, 'YLim');
            ylim([0 yy(2)])
            set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
            if i==1
                mes = sprintf('psth \n (pref - null ch)');
                ylabel(mes)
            else
                set(gca, 'XTick', [])
            end
            
            % PK
            subplot(7,repeat,i+6*repeat)
            norm = mean([yout.amp_h, yout.amp_l]);
            plot(1:7, yout.amp_h/norm, '-','color',y)
            hold on;
            plot(1:7, yout.amp_l/norm, '-','color',g)
            xlim([1 7])
            yy = get(gca, 'YLim');
            ylim([0 yy(2)])
            set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
            if i==1
                ylabel('PK')
            else
                set(gca, 'XTick', [])
            end
            
            disp('successfully run...')
            disp('------------------------------')
%         catch
%             disp('error occured...skipped.')
%         end

    end
    
    switch k
        case 1
            figname = 'contrast kernel: variable, stimlus kernel: 0.015';
        case 2
            figname = 'contrast kernel: 0.015, stimulus kernel: variable';
    end
    set(gcf, 'Name', figname, 'NumberTitle', 'off')
end