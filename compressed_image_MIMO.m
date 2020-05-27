clc; clear; close all;

img = double(rgb2gray(imread('lena_std.tif')));

[U,S,V]=svd(img);

[img_row,img_column] = size(img);

rank = 150;

u1 = U(:,1:rank); s1 = diag(S(1:rank,1:rank)); v1 = V(:,1:rank);

data = [u1;v1;s1'];

bit = [0 1];

point1 = nrSymbolModulate(bit','BPSK'); s0 = point1(1); s1 = point1(2);

cons(1,1) = s0; cons(2,1) = s1;

point = [s0  s0 ; s0 s1 ; s1 s0 ; s1 s1];

[a, b] = size(data);

rounded_data = round(800*data(1:a-1,:))+200;

m_bi1 = de2bi(rounded_data); %data를 2진수로 encoding

m_bi2 = de2bi(round(data(a,:)));

[r1,c1] = size(m_bi1); [r2,c2] = size(m_bi2);

m_bi_vec1 = reshape(m_bi1,1,r1*c1); m_bi_vec2 = reshape(m_bi2,1,r2*c2); 

m_bi_vec = [m_bi_vec1,m_bi_vec2];

n = 1 : size(m_bi_vec,2);

symbol = nrSymbolModulate(m_bi_vec','BPSK');

sym = reshape(symbol,length(m_bi_vec)/2,2);

mimochannel = comm.MIMOChannel('PathGainsOutputPort', true ); % default 값 Tx 2개 Rx 2개

[recieve_signal,H] = step(mimochannel,sym);

re = reshape(recieve_signal,length(recieve_signal)*2,1); % 노이즈 더할려고 폄

output_ZF = zeros(length(recieve_signal),2);

EbNo = -10 : 2 : 20;

CSI = zeros(2,2);

ber_ZF = zeros(1,length(EbNo));

ber_ML = zeros(1,length(EbNo));

ber_MMSE = zeros(1,length(EbNo));

k = zeros(4,2); like = zeros(1,4); 

output_ML = zeros(length(recieve_signal),2);

output_MMSE = zeros(length(recieve_signal),2);

for j = 1 : length(EbNo)

   
    re_awgn = awgn(re,EbNo(j),'measured');
    
    recieve = reshape(re_awgn,length(re_awgn)/2,2);
    
    for i = 1 : length(recieve)  %Zero Forcing 
    
        %CSI = squeeze(H(i,1,:,:)); 연산 너무 느림 
        
        CSI(1,1) = H(i,1,1,1); CSI(1,2) = H(i,1,1,2); CSI(2,1) = H(i,1,2,1); CSI(2,2) = H(i,1,2,2);
        
        h = recieve(i,:) * pinv(CSI);
        
        output_ZF(i,1) = h(1,1);
        
        output_ZF(i,2) = h(1,2);
        
        for l = 1 :4  %ML
            k(l,:) = point(l,:)*CSI;
        end
          
        for l = 1 :4
            like(1,l) = norm(recieve(i,:)-k(l,:));
        end
        [q,lo] = min(like);
        
        output_ML(i,:) = point(lo,:);
        
        % MMSE
        No =sqrt(2)/ (10^(EbNo(j)/10));
        
        W = pinv(CSI'*CSI + No*eye(2))*CSI';
        
        mm = recieve(i,:) * W;
        
        output_MMSE(i,1) = mm(1,1);
        
        output_MMSE(i,2) = mm(1,2);
        
    end
    %output_SD = sphDec(recieve,squeeze(H));
        
    ZF = reshape(output_ZF,length(output_ZF)*2,1);
        
    ML = reshape(output_ML,length(output_ML)*2,1);
    
    MMSE = reshape(output_MMSE,length(output_MMSE)*2,1);
        
    demod_ZF = nrSymbolDemodulate(ZF,'BPSK');
        
    demod_ML = nrSymbolDemodulate(ML,'BPSK');
    
    demod_MMSE = nrSymbolDemodulate(MMSE,'BPSK');
        
    for l = 1 : length(demod_ZF)
        if demod_ZF(l)<0
            demod_ZF(l) = 1;
        else
            demod_ZF(l) = 0;
        end
    end
        
    for l = 1 : length(demod_ML)
        if demod_ML(l)<0
            demod_ML(l) = 1;
        else
            demod_ML(l) = 0;
        end
    end
    
    for l = 1 : length(demod_MMSE)
       if demod_MMSE(l)<0
            demod_MMSE(l) = 1;
       else
            demod_MMSE(l) = 0;
       end
    end
    
    k_ZF = abs(demod_ZF - m_bi_vec');

    ber_ZF(j) = sum(k_ZF)/length(demod_ZF);
    
    k_ML = abs(demod_ML - m_bi_vec');

    ber_ML(j) = sum(k_ML)/length(demod_ML);
    
    k_MMSE = abs(demod_MMSE - m_bi_vec');

    ber_MMSE(j) = sum(k_MMSE)/length(demod_MMSE);
    
    semilogy(EbNo , ber_ZF,'ro')
    
    hold on;
    
    semilogy(EbNo , ber_ML,'b*')
    
    semilogy(EbNo , ber_MMSE,'g+')
    
    legend('Zero Forcing','Maximum Likelihood','L-MMSE')
    
    grid on; 
    
    xlim([-10, 20]);
    
    ylim([1e-5, 1]); 
    
    xlabel('Eb/No (dB)');
    
    ylabel('BER');
    
    legend('Zero Forcing','Maximum Likelihood','L-MMSE')
    
    drawnow;
    
    
        
end

hold off;

reshape_MMSE_bi = reshape(demod_MMSE(1:length(m_bi_vec1)),r1,c1);

singular_MMSE_bi = reshape(demod_MMSE(length(m_bi_vec1)+1:length(demod_MMSE)),r2,c2);

reshape_MMSE = bi2de(reshape_MMSE_bi); singular_MMSE = bi2de(singular_MMSE_bi);

re_MMSE = [(reshape(reshape_MMSE ,a-1,b)-200);singular_MMSE']; 

u = re_MMSE(1:img_row,:); 

v = re_MMSE(img_row+1:img_row+img_column,:); 

s = diag(re_MMSE(img_row+img_column+1,:));

figure;

imagesc(u*s*v')

colormap('gray')







