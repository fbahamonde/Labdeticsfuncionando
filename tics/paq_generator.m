function y = paq_generator(imagen,bpS,d_bits)
    r = length(imagen(:,1));
    c = length(imagen(1,:));
    Ns = ceil(c*r/d_bits);
    l_gen = zeros(1,Ns);
    data_vector = reshape(imagen,1,c*r);
    bits_extras = (Ns*d_bits - c*r);
    if bits_extras==0 
         for i = 1:Ns
             data = data_vector(1+(i-1)*d_bits:i*d_bits);
             b_sum = sum(data)
             p_bit = mod((b_sum),2)
             p_data = horzcat(data,p_bit);
             paq_list{i} = p_data;
        end
    else
         for i = 1:Ns-1
             data = data_vector(1+(i-1)*d_bits:i*d_bits);
             b_sum = sum(data);
             p_bit = mod((b_sum),2);
             p_data = horzcat(data,p_bit);
             paq_list{i} = p_data;
         end   
         paq_extra = zeros(1,bpS-1);
         n_extrab = d_bits - bits_extras; %Número de data bits en el paquete extra
         paq_extra(1:n_extrab) = data_vector(end+1-n_extrab:end);
         b_sum = sum(paq_extra);
         p_bit = mod((b_sum),2);
         p_data = horzcat(paq_extra,p_bit);
         paq_list{Ns} = p_data;
    end     
    y = paq_list;
end