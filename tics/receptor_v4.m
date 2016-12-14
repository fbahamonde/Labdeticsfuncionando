clear
close all
clc

%% Simulation Display 
plot_switch = 'off';
sound_switch = 'on';
place_fig = 16;
place_fig_v = 478; % 1558;
fig_width = 480; %640;
fig_height = 360; %480;
max_width = 1920; %3840;
max_height = 964; %2044;

%% Initialization
Total = 17;
counter = 1;
bpS=16;
d_bits=15;
N_paq = 2; %Estado inicial
paq_actual = 1;
done = 0;
fs = 40960;         % sampling frequency
N_frame = 0;        % samples in the frame
t_buffer = 1;       % buffer time
t_rec = 5;   % recording time
t_sync = 0.1; % duration
t_ch = 0:1/fs:t_sync-1/fs;
f0_ch = fs/2^4;
f1_ch = fs/2^3;
z_sync = 16;
sync = [chirp(t_ch(1:end-z_sync),f0_ch,t_sync-z_sync/fs,f1_ch,'linear',90),zeros(1,z_sync)];
t_s = 0.25;                                     % symbol period
f_start = 2688;                                 % start frequency
f_end = 17600; %16064; % 11904;                                  % end frequency
f_step = 800;
f_end = f_start + 17*f_step;
t = 0:1/fs:t_s-1/fs;                            % time vector (symbol)
Nt = length(t);                                 % total samples in symbol
f = 0:fs/Nt:fs*(1-1/Nt);                        % frequency vector
Nf = length(f);                                 % total frequency samples in the symbol
fc = f_start:f_step:f_end;                      % carrier frequencies
Nfc = length(fc);                               % amount of carrier frequencies
fc_index = fc*t_s + 1;                          % frequency vector indexes of the carrier frequencies

disp('Presione una tecla para continuar ')
pause()
while done == 0 
contador2=1;
fprintf('Receiving ... ')
h_rx = audiorecorder(fs,16,1);
recordblocking(h_rx,t_rec);
rx_frame = getaudiodata(h_rx, 'double')';

% figure placement reset
place_fig = 16;
place_fig_v = place_fig_v - 92 - fig_height;

%% Leer imagen
pic = imread('cara.jpg');
pic = im2bw(pic);
paquetes = paq_generator(pic,bpS,15);
[c,r] = size(pic);
bpS = Nfc - 2;  
Ns = ceil(c*r/d_bits);
data = zeros(1,Nt*Ns);



%% Synchronization
tic
fprintf('Finding signal starting point ... ')
find_max = conv(rx_frame(1:fs), fliplr(sync));
[~, n_start] = max(find_max);

switch plot_switch
    case 'on'
        figure('Name','Start Point','NumberTitle','off','Color',[0.75 0.75 0.75],'Position',[place_fig place_fig_v fig_width fig_height])
        plot(rx_frame(1:fs))
        hold on
        plot(find_max/max(find_max),'c')
        plot(n_start+1,rx_frame(n_start+1),'or')
        place_fig = place_fig + 16 + fig_width;
end
elapsed_time = toc;
fprintf('done. [%f seconds]\n',elapsed_time)


%% Recover Signal
tic
fprintf('Computing transfer function ... ')
rx_ch = 2*rx_frame(n_start+1:n_start+Nt);
Channel = fft(rx_ch);


%% Transfer Function
TFunc = abs(Channel(fc_index))*2*Nfc/Nf;
Fix = TFunc.^-1;
Fix_TFunc = [zeros(1,f_start*Nt/fs), upsample(Fix,f_step*Nt/fs), zeros(1,Nf-f_start*Nt/fs-Nfc*f_step*Nt/fs)];
Channel2 = Channel.*Fix_TFunc;

switch plot_switch
    case 'on'
        figure('Name','Received Channel Signal','NumberTitle','off','Color',[0.75 0.75 0.75],'Position',[place_fig place_fig_v fig_width fig_height])
        subplot(3,1,1)
        plot(t, rx_ch)
        subplot(3,1,2)
        plot(f, abs(Channel2)*2*Nfc/Nf)
        subplot(3,1,3)
        plot(f, angle(Channel2).*abs(Channel2)*2*Nfc/Nf)
        place_fig = place_fig + 16 + fig_width;
end
elapsed_time = toc;
fprintf('done. [%f seconds]\n',elapsed_time)


%% Recover Pilot
tic
fprintf('Recovering pilot signal ... ')
n_start = n_start + Nt;
rx_pilot = 2*rx_frame(n_start+1:n_start+Nt);
Pilot = fft(rx_pilot).*Fix_TFunc;

switch plot_switch
    case 'on'
        figure('Name','Received Pilot Signal','NumberTitle','off','Color',[0.75 0.75 0.75],'Position',[place_fig place_fig_v fig_width fig_height])
        subplot(3,1,1)
        plot(t, rx_pilot)
        subplot(3,1,2)
        plot(f, abs(Pilot)*4/Nf)
        subplot(3,1,3)
        plot(f, angle(Pilot).*abs(Pilot)*4/Nf)
        place_fig = place_fig + 16 + fig_width;
end
elapsed_time = toc;
fprintf('done. [%f seconds]\n',elapsed_time)


%% Recover Data
tic
fprintf('Recovering data ... ')
threshold = 0.1;
n_start = n_start + Nt;
rx_data = zeros(N_paq, bpS);
%figure('Name','Data Retrieval','NumberTitle','off','Color',[0.75 0.75 0.75],'Position',[place_fig place_fig_v fig_width fig_height])
place_fig = place_fig + 16 + fig_width;
for n = 1:N_paq
    rx_data_signal = 2*rx_frame(n_start+1:n_start+Nt);
    Data = fft(rx_data_signal).*Fix_TFunc;
    fc_data = abs(Data(fc_index));
    fc_data = fc_data/max(fc_data);
    pause(0.01)
    rx_fc_data = fc_data > threshold;
    rx_data(n,:) = rx_fc_data(2:end-1);
    n_start = n_start + Nt;
end
N_paq;
rx_data;

switch plot_switch
end
elapsed_time = toc;
fprintf('done. [%f seconds]\n',elapsed_time)
p_rx{counter}=rx_data;
paq_actual = paq_actual + N_paq;

%% retransmisión y bit de paridad
prx1=p_rx{counter};
counter = counter + 1;
e_data = zeros(1,N_paq+1);
e_data(1) = N_paq;

for i = 1:N_paq
    p = prx1 (i,:);
    suma = sum(p(1:15));
    xor = mod(suma,2);
    if xor ~= p(16)
        e_data(i+1) = 1;
    end        
end
edatas{contador2}=e_data;
edat=edatas{contador2};
err_perc = sum(e_data(2:end))/edat(1);
contador2=contador2+1;
umbral = 0.5;
if err_perc > umbral
   N_paq = ceil(N_paq/2);
   counter=counter-1;
   paq_actual=paq_actual-edat(1);
   disp('Pidiendo retransmisión')
else
    N_paq = N_paq + 1;
    disp('Se pasa al siguiente paquete')
end
if paq_actual + N_paq > 34
    N_paq = 34 + 1 - paq_actual;
end    
if paq_actual > 34
    done = 1;
end    

%% Conexión cliente
t = tcpip('192.168.43.130',4013);
fopen(t);
pause(0.2);
fwrite(t,e_data);
fclose(t);
pause(3-0.1)

end


%% recuperando imagen
unionn=[];
for i=1:length(p_rx);
    unionn=vertcat(unionn,p_rx{i});
end
rr=size(unionn);
union2=zeros(34,15);

for i=1:34
    for j=1:15
        union2(i,j)=unionn(i,j);
    end
end
union3=reshape(union2',1,34*15);
for i=1:21*24
    union4(i)=union3(i);
end
union5=reshape(union4,21,24);
recibida=union5;
true_img=rgb2gray(imread('cara.jpg'));
true_img=true_img>128;
figure (1)
subplot(1,2,1)
imshow(recibida);
xlabel('Recibida')
subplot(1,2,2)
imshow(true_img);
xlabel('Original')