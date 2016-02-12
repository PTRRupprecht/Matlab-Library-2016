function LFP_analysis
% ECS: external command sensitivity [mV/V]
% DAQfreq: data acquisition frequency [Hz]
% FSP 2008/04/02

global ECS
ECS = 1;                %scaling factor
global DAQfreq
DAQfreq = 10000;
global rootpath
rootpath = 'C:\data\sbund\';
global G_FileNames
G_FileNames = {''};
%G_FileNames = {'s0003080905'; 's0001080907'; 's0002080907'; 's0001080909'; 's0003080909'; 's0001080910'};
G_FileNames = {'s100410'};
%G_FileNames = {'s100121'};
global G_StartEnd
G_StartEnd = [17, 54]; % [7, 43];
global G_Bandwidth
G_Bandwidth = [5 40];
global G_TimeWindow
G_TimeWindow = [1,4, 4, 10];
global G_FreqBand
G_FreqBand = [5, 40];
global LFP_trace
LFP_trace = 1;
%--------------------------------------------------------------------------          

maxLength = max(G_StartEnd(:,2) - G_StartEnd(:,1)) + 1;
PowerMat = zeros(maxLength, size(G_FileNames,1));

for ii = 1 : size(G_FileNames, 1)
    
    size(G_FileNames, 2);
    
    tmpLength = G_StartEnd(ii,2) - G_StartEnd(ii,1) + 1;
    tmpPower = zeros(tmpLength, 1);
    xScale_time = zeros(tmpLength, 1);

    L_Basename = [rootpath, G_FileNames{ii}, '\sb0001\sb0001AAAA'];
    
    jjcounter = 0;
    
    for jj = G_StartEnd(ii, 1) : G_StartEnd(ii, 2)
    
        jjcounter = jjcounter + 1;
        
        if (jj > 999)

            L_FileName = [L_Basename, num2str(jj), '.xsg'];

            elseif (jj > 99)
            
                L_FileName = [L_Basename, '0', num2str(jj), '.xsg'];
                    
            elseif (jj > 9)
                    
                L_FileName = [L_Basename, '00', num2str(jj), '.xsg'];
                    
            else
                    
                L_FileName = [L_Basename, '000', num2str(jj), '.xsg'];
                    
        end
                
        L_Data = load(L_FileName, '-mat');
        if (LFP_trace == 1)
            L_trace = L_Data.data.ephys.trace_1;
        elseif (LFP_trace == 2)
            L_trace = L_Data.data.ephys.trace_2;
        end
        L_time = L_Data.header.xsgFileCreationTimestamp(13:20);
        xScale_time(jjcounter) = str2num(L_time(1:2))*60 + str2num(L_time(4:5)) + str2num(L_time(7:8))/60;
        
        
        if (jj == G_StartEnd(ii, 2))
            
            xdata = [1 : size(L_trace)] / DAQfreq; %in sec
         
            new_data = timeseries(L_trace,xdata);
            idealfilter_data = idealfilter(new_data,G_Bandwidth,'pass');
            figure
            plot(new_data,'b-'), hold on
            plot(idealfilter_data,'r-')
            hold off
            
        end
        
        L_trace_baseline = L_trace(G_TimeWindow(1) * DAQfreq : G_TimeWindow(2) * DAQfreq);
        L_trace_response = L_trace(G_TimeWindow(3) * DAQfreq : G_TimeWindow(4) * DAQfreq);
        freq_norm = G_FreqBand(2) - G_FreqBand(1);
        
        [Pxx,f] = pwelch(L_trace_baseline,length(L_trace_baseline),[],length(L_trace_baseline),DAQfreq);
        fmin_pos = min(find(f > G_FreqBand(1)));
        fmax_pos = max(find(f < G_FreqBand(2)));
        tmpPower(jjcounter,1) = sum(Pxx(fmin_pos:fmax_pos))/freq_norm;
        
        [Pxx,f] = pwelch(L_trace_response,length(L_trace_response),[],length(L_trace_response),DAQfreq);
        fmin_pos = min(find(f > G_FreqBand(1)));
        fmax_pos = max(find(f < G_FreqBand(2)));
        tmpPower(jjcounter,2) = sum(Pxx(fmin_pos:fmax_pos))/freq_norm;
        
    end
    
    xScale_Trial = [G_StartEnd(ii, 1) : G_StartEnd(ii, 2)];
    figure
    plot(xScale_Trial,tmpPower(:,1),'ko-');
	hold on
    plot(xScale_Trial,tmpPower(:,2), 'ro-');
    hold off

    xScale_time = xScale_time - xScale_time(1);
    figure
    plot(xScale_time,tmpPower(:,1),'ko-');
	hold on
    plot(xScale_time,tmpPower(:,2), 'ro-');
    hold off

end



