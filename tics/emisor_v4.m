clear
close all
clc
done = 0; %No he terminado de enviar la imagen.
link = 0;
ready = 1;
d_bits = 15;

%CARA
pic = imread('cara.jpg');
pic = im2bw(pic); %binarización
[c,r] = size(pic);
bpS = 16;  % bits per symbol
Ns = ceil(c*r/d_bits);


%Vuelta al codigo original
paquetes = paq_generator(pic,bpS,15);
N_paq = 2; %Se inicia con 2 paquetes por envío.
paq_actual = 1; 
k=1; %Indice de acks recibidos
pause()

while done == 0    
    %pause()
    while ready ==1
        k
        %paq_actual = paq_actual+2        
        %% Simulation Display
        
        s = SimDisp('off','on',1); % Argumentos: sound switch, plot switch y modo (1 o 2, determina valor de parametros)
        plot_switch = s(1);
        sound_switch = s(2);
        
        %% Initialization
        fs = 40960;         % sampling frequency
        N_frame = 0;        % samples in the frame
        t_buffer = 1;       % buffer time
        t_rec = t_buffer;   % recording time
        
        
        %% Synchronization
        tic
        %fprintf('Generating sync pulse ... ')
        t_sync = 0.1; % duration
        t_ch = 0:1/fs:t_sync-1/fs;
        f0_ch = fs/2^4;
        f1_ch = fs/2^3;
        z_sync = 16;
        sync = [chirp(t_ch(1:end-z_sync),f0_ch,t_sync-z_sync/fs,f1_ch,'linear',90),zeros(1,z_sync)];
        switch plot_switch
            case 'on'
                figure('Name','Synchronization','NumberTitle','off','Color',[0.75 0.75 0.75],'Position',[place_fig place_fig_v fig_width fig_height])
                plot(t_ch,sync,'.-')
                hold on
                plot(t_ch(end-z_sync),sync(end-z_sync),'or')
                place_fig = place_fig + 16 + fig_width;
                %         axis([0.09,0.11,-1,1])
        end
        N_sync = length(sync);
        N_frame = N_frame + N_sync;
        t_rec = t_rec + t_sync;
        elapsed_time = toc;
        %fprintf('done. [%f seconds]\n',elapsed_time)
        
        
        %% Carrier Frequencies
        tic
        %fprintf('Generating channel characterization signal ... ')
        t_s = 0.25;                                     % symbol period
        f_start = 2688;                                 % start frequency
        f_end = 17600; %16064; % 11904;                                  % end frequency
        %frec hack
        %f_step=ceil((f_end-f_start)/201);
        %Homero normal y traspuesto
        %f_step=64;  %
        %f_end=f_start+f_step*201; % 201 o 233 , trasp o normal
        %f_end=f_start+f_step*201;
        %CARA
        f_step=800;
        f_end = f_start +  f_step*17;
        %fin hack
        %f_step = 64;                                    % frequency step (separation between carrier frequencies)
        t = 0:1/fs:t_s-1/fs;                            % time vector (symbol)
        Nt = length(t);                                 % total samples in symbol
        f = 0:fs/Nt:fs*(1-1/Nt);                        % frequency vector
        Nf = length(f);                                 % total frequency samples in the symbol
        fc = f_start:f_step:f_end;                      % carrier frequencies
        Nfc = length(fc);                               % amount of carrier frequencies
        fc_index = ceil(fc*t_s) + 1;                          % frequency vector indexes of the carrier frequencies
        
        
        %% Channel Characterization Signal
        ch_signal = zeros(1, Nt);
        for fi = fc
            ch_signal = ch_signal + sin(2*pi*fi*t);
        end
        ch_signal = 1/Nfc*ch_signal;                    % normalize signal
        CH_Signal = fft(ch_signal);                     % frequency content of ch_signal
        
        N_frame = N_frame + Nt;
        t_rec = t_rec + t_s;
        
        switch plot_switch
            case 'on'
                figure('Name','Sent Channel Signal','NumberTitle','off','Color',[0.75 0.75 0.75],'Position',[place_fig place_fig_v fig_width fig_height])
                subplot(3,1,1)
                plot(t,ch_signal)
                subplot(3,1,2)
                plot(f, abs(CH_Signal)*2*Nfc/Nf)
                subplot(3,1,3)
                plot(f, angle(CH_Signal).*abs(CH_Signal)*2*Nfc/Nf)
                place_fig = place_fig + 16 + fig_width;
        end
        elapsed_time = toc;
        %printf('done. [%f seconds]\n',elapsed_time)
        
        
        %% Pilot Signal
        tic
        %fprintf('Generating pilot signal ... ')
        pilot_signal = 0.5*(sin(2*pi*fc(1)*t) + sin(2*pi*fc(end)*t));
        Pilot_Signal = fft(pilot_signal);                     % frequency content of ch_signal
        
        N_frame = N_frame + Nt;
        t_rec = t_rec + t_s;
        
        switch plot_switch
            case 'on'
                figure('Name','Sent Pilot Signal','NumberTitle','off','Color',[0.75 0.75 0.75],'Position',[place_fig place_fig_v fig_width fig_height])
                subplot(3,1,1)
                plot(t,pilot_signal)
                subplot(3,1,2)
                plot(f, abs(Pilot_Signal)*4/Nf)
                subplot(3,1,3)
                plot(f, angle(Pilot_Signal).*abs(CH_Signal)*4/Nf)
                place_fig = place_fig + 16 + fig_width;
        end
        elapsed_time = toc;
        %fprintf('done. [%f seconds]\n',elapsed_time)
        
        
        
        %% Data
        tic
        fprintf('Generating data signal ... ')
        
        image_bin = zeros(N_paq,bpS);
        for i = 1:N_paq           
            image_bin(i,:) = paquetes{paq_actual-1+i};
        end
        %image_bin = reshape(pic, Ns, bpS);
        data_matrix = [ones(N_paq,1), image_bin, ones(N_paq,1)];
        
        data = zeros(1,Nt*Ns);
        
        for n = 1:N_paq
            Nfn = sum(data_matrix(n,:));
            for fi = 1:Nfc
                if data_matrix(n,fi)
                    data(1+Nt*(n-1):Nt*n) = data(1+Nt*(n-1):Nt*n) + sin(2*pi*fc(fi)*t);
                end
            end
            data(1+Nt*(n-1):Nt*n) = 1/Nfn*data(1+Nt*(n-1):Nt*n);
        end
        
        t_data = N_paq*t_s;
        t_rec = t_rec + t_data;
        N_data = Nt*Ns;
        N_frame = N_frame + N_data;
        elapsed_time = toc;
        fprintf('done. [%f seconds]\n',elapsed_time)
        
        
        %% Build whole signal
        tic
        fprintf('Building whole signal ... ')
        tx_frame = [sync,ch_signal,pilot_signal,data];
        t_frame = 0:1/fs:(N_frame - 1)/fs;
        
        switch plot_switch
            case 'on'
                figure('Name','Frame','NumberTitle','off','Color',[0.75 0.75 0.75],'Position',[place_fig place_fig_v fig_width fig_height])
                plot(t_frame, tx_frame)
                %         place_fig = place_fig + 16 + fig_width;
        end
        elapsed_time = toc;
        fprintf('done. [%f seconds]\n',elapsed_time)
        
        
        %% Transmit using Sound
        tic
        fprintf('Transmitting ... ')
        % Transmitting Frame
        %pause()
        h_tx = audioplayer(tx_frame,fs);
        play(h_tx);
        paq_actual = paq_actual + N_paq;
        %image_bin
        k
        N_paq        
        %SERVER TCP PARA RECIBIR CONFIRMACIONES
        disp ('Receiver started');
        t=tcpip('0.0.0.0', 4013,'NetworkRole','server','timeout',3,'inputbuffersize',504);
        
        % Wait for connection
        disp('Waiting for connection');
        fopen(t);
        pause(0.2);
        disp('Connection OK');
        tic
        f_data{k}=fread(t);
        k=k+1;        
        fclose(t);
        toc
        ready = 0;
    end
    %%Confirmación de recepción y errores.
    rx_paq = f_data{k-1};
    n_prev = rx_paq(1);  %numero de paquetes enviados previamente
    err_perc = sum(rx_paq(2:end))/n_prev;
    umbral = 0.5;
    if err_perc > umbral
        re = 1; %Existe retransmisión
        N_paq = ceil(N_paq/2); %ceil evita enviar 0 paquetes
        paq_actual = paq_actual - n_prev;
        disp('Retransmitiendo...')
    else 
        N_paq = N_paq + 1;
    end
    %Envío paquetes de sobra?
    if paq_actual + N_paq > 34
        %Envío los que me faltan
        N_paq = 34 + 1 - paq_actual;
    end    
    %Envié todo?
    if paq_actual >34
        done = 1;
    end
    ready = 1;
end