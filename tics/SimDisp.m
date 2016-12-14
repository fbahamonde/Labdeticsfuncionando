function y = SimDisp(p_switch,s_switch,modo)
    plot_switch = p_switch ;
    sound_switch = s_switch ;
    if modo==1        
    place_fig = 16;
    place_fig_v = 478; % 1558;
    fig_width = 480; %640;
    fig_height = 360; %480;
    max_width = 1920; %3840;
    max_height = 964; %2044;    
    else
    place_fig = 16/2;
    place_fig_v = 1558/2;
    fig_width = 640;
    fig_height = 480;
    max_width = 3840/2;
    max_height = 2044/2;
    end
    y = [p_switch s_switch];
end

