close all
clear

% Paramètres de l'étude
Rb = 3000;      % débit binaire
Fe = 24000;     % fréquence d'échantillonnage
N = 20000;      % nombre de message_binaire à transmettre
M = 4;          % ordre de modulation
Rs = Rb/log2(M);  % débit symbole
Ns = sqrt(Fe /Rs); % taux d'échantillonnage
alpha = 0.35; % roll off
L = 8;

% Génération du message binaire aléatoire
n = 20000; % Nombre de message_binaire qu'on veut transmettre
message_binaire = randi([0, 1], 1, n);

%% Codage canal

% Introduction du code convolutif 
trellis = poly2trellis(7,[171 133]);  % Initialisation de la structure de code
bits_codes = convenc(message_binaire, trellis);  % Codage des message_binaire en utilisant la fonction "convenc"
N_bits_codes = length(bits_codes);    % Nombre de message_binaire codés


% Poinçonnage
matrice_p = [1 1 0 1];          % Matrice de poinçonnage
bits_poinconnes = convenc(message_binaire, trellis, matrice_p);
N_poinconnes = length(bits_poinconnes);


%% Sans poinçonnage

% Mapping
symboles_codes = qammod(bits_codes.',M,'InputType','bit').';

% Suréchantillonnage
s_codes = kron(symboles_codes, [1 zeros(1, Ns-1)]);
h = rcosdesign(alpha, L, Ns, 'sqrt');
signal_e_codes = filter(h, 1, [s_codes zeros(1, length(h)-1)]);


% Canal AWGN passe-bas équivalent pour le cas sans poinçonnage
R = 0.5;          % Rendement du code sans poinçonnage 
Eb_N0_dB2 = -10:12;
%Eb_N0_dB = 4:7;  % valeurs SNR (Eb/No) en (dB) allant de 4 à 7 dB pour comparer les TEB

Eb_N02 = 10.^(Eb_N0_dB2/10);
P_signal_e_codes = mean(abs(signal_e_codes).^2); % Puissance du signal émis
sigma_ne_codes = sqrt((R*P_signal_e_codes * Ns)./(2*log2(M)*Eb_N02));
dim1 = length(Eb_N02);
dim2 = length(signal_e_codes);

n_I_codes = randn(dim1, dim2) .* sigma_ne_codes.';
n_Q_codes = randn(dim1, dim2) .* sigma_ne_codes.';
n_e_codes = n_I_codes + 1i * n_Q_codes;

r_codes = repmat(signal_e_codes, dim1, 1) + n_e_codes;


% Réception
h_reception = h;
z_codes = filter(h_reception,1, r_codes,[],2);

% Echantillonnage
n0 = length(h);
z_codes = z_codes(:, n0: Ns :end);


% Décision
d_codes = sign(real(z_codes)) + 1i * sign(imag(z_codes));

bits_recus_codes = zeros(dim1, N_bits_codes);
bits_recus_soft = zeros(dim1, N_bits_codes);

for i = 1:dim1
    bits_recus_codes(i,:) = (qamdemod(d_codes(i,:).', M, 'OutputType', 'bit')).';
    bits_recus_soft(i, :) = (qamdemod(z_codes(i,:).', M, 'OutputType', 'llr')).';
end
clear i;

% Décodage canal
tblen = 5*6;
bits_decodes_hard = zeros(dim1, N);
bits_decodes_soft = zeros(dim1, N);
for i=1:dim1
    bits_decodes_hard(i,:) = vitdec(bits_recus_codes(i,:),trellis,tblen,'trunc','hard');
    bits_decodes_soft(i,:) = vitdec(bits_recus_soft(i,:),trellis,tblen,'trunc','unquant');
end

% Calcul du taux d'erreur binaire
TEB_soft = mean(abs(message_binaire-bits_decodes_soft),2).';
TEB_hard = mean(abs(message_binaire-bits_decodes_hard),2).';

figure()
semilogy(Eb_N0_dB2, TEB_hard);
grid on;
hold on
semilogy(Eb_N0_dB2, TEB_soft);
grid on;
hold on
semilogy(Eb_N0_dB2, qfunc(sqrt(2*Eb_N02)));
legend("Hard R=1/2","Soft R=1/2", 'Théorique')
xlabel('E_b/N_0 (dB)')
ylabel('TEB')
title('Comparaison des TEBs avec codage canal')


%% Avec poinçonnage

% Mapping  
symboles_poinconne = qammod(bits_poinconnes.',M,'InputType','bit').';

% Suréchantillonnage
s_ap = kron(symboles_poinconne, [1 zeros(1, Ns-1)]);
signal_e_poinconne = filter(h, 1, [s_ap zeros(1, length(h)-1)]);

% Canal AWGN passe-bas équivalent pour le cas avec poinçonnage
R_poinconnage = 2/3;       % Rendement du code avec poinçonnage
P_signal_e_poinconne = mean(abs(signal_e_poinconne).^2); % Puissance du signal émis
sigma_ne_poinconne = sqrt((R_poinconnage*P_signal_e_poinconne * Ns)./(2*log2(M)*Eb_N02));
dim1 = length(Eb_N02);
dim2 = length(signal_e_poinconne);
n_I_poinconne = randn(dim1, dim2) .* sigma_ne_poinconne.';
n_Q_poinconne = randn(dim1, dim2) .* sigma_ne_poinconne.';
n_e_poinconne = n_I_poinconne + 1i * n_Q_poinconne;


r_poinconne = repmat(signal_e_poinconne, dim1, 1) + n_e_poinconne;

% Réception
z_poinconne = filter(h_reception,1, r_poinconne,[],2);

% Echantillonnage
z_poinconne = z_poinconne(:, n0: Ns :end);

% Démodulation
d_poinconne = sign(real(z_poinconne)) + 1i * sign(imag(z_poinconne));

% Décision
bits_recus_poinconne = zeros(dim1, N_poinconnes);
for i = 1:dim1
    bits_recus_poinconne(i,:) = (qamdemod(d_poinconne(i,:).', M, 'OutputType', 'bit')).';
end
clear i;

% Décodage canal
tblen = 5*7;
bits_decoded_ap = zeros(dim1, N);
for i=1:dim1
    bits_decoded_ap(i,:) = vitdec(bits_recus_poinconne(i,:),trellis,tblen,'trunc','hard', matrice_p);
end

% Calcul du taux d'erreur binaire
TEB_poinc = zeros(1,length(Eb_N02));
for i=1:length(Eb_N02)
    TEB_poinc(i) = mean(abs(message_binaire-bits_decoded_ap(i,:)));
end
clear i;
 
% Comparaison entre TEB avec poinçonnage et TEB sans poinçonnage

figure()
semilogy(Eb_N0_dB2, TEB_hard,'--');
grid on;
hold on
semilogy(Eb_N0_dB2, TEB_poinc);
grid on;
legend("sans poinçonnage R=1/2","avec poinçonnage R=2/3")
xlabel('E_b/N_0 (dB)')
ylabel('TEB')
title('Comparaison entre TEB avec poinçonnage et TEB sans poinçonnage')

%% Ajout de l'entrelaceur et du code de Reed Solomon
% Paramètres de l'entrelaceur
N = 204;         % La longueur totale des codes RS(N)
K = 188;         % La longueur des données originales avant codage RS(K)

% Génération de bits
bits_K = randi([0 1], 1, K*8*100);

% Entrelaceur
encoder = comm.RSEncoder(N, K, 'BitInput', true); % encodeurr de la classe RSEncoder pour effectuer le codage RS(N,K)
decoder = comm.RSDecoder(N, K, 'BitInput', true); % decodeur de la classe RSEncoder pour effectuer le decodage RS(N,K)
bitscoderRS = step(encoder, bits_K.');
vect_perm = randperm(length(bitscoderRS)); % vect_perm de permutations aléatoires de la longueur de bitscoderRS
intr = reshape(bitscoderRS(vect_perm), [], 1);


% Codage canal
bits_K_codes = convenc(bits_K, trellis, matrice_p);
bits_K_codes_RS = convenc(intr.', trellis, matrice_p);

N_K = length(bits_K_codes);
N_K_RS = length(bits_K_codes_RS);

% Modulation et filtrage de mise en forme
symboles_K_codes = (qammod(bits_K_codes.', M, 'InputType', 'bit')).';
symboles_K_codes_RS = (qammod(bits_K_codes_RS.', M, 'InputType', 'bit')).';


s_K_codes = kron(symboles_K_codes, [1 zeros(1, Ns-1)]);
s_K_codes_RS = kron(symboles_K_codes_RS, [1 zeros(1, Ns-1)]);
signal_e_K_codes = filter(h, 1, [s_K_codes zeros(1, length(h)-1)]);   % Signal émis sans codage RS
signal_e_K_codes_RS = filter(h, 1, [s_K_codes_RS zeros(1, length(h)-1)]); % Signal émis avec codage RS

% Canal AWGN passe-bas équivalent
R_K = 2/3; % Rendement du code 
P_signal_e_K_codes = mean(abs(signal_e_K_codes).^2); % Puissance du signal émis
sigma_ne_K_codes = sqrt((R_K*P_signal_e_poinconne * Ns)./(2*log2(M)*Eb_N02));
dim1 = length(Eb_N02);
dim2 = length(signal_e_K_codes);
n_I_K_codes = randn(dim1, dim2) .* sigma_ne_K_codes.';
n_Q_K_codes = randn(dim1, dim2) .* sigma_ne_K_codes.';
n_e_K_codes = n_I_K_codes + 1i* n_Q_K_codes;

r_K_codes = repmat(signal_e_K_codes, dim1, 1) + n_e_K_codes;


% Canal AWGN passe-bas équivalent
R_K_RS = 2/3; % Rendement du code
P_signal_e_K_codes_RS = mean(abs(signal_e_K_codes_RS).^2); % Puissance du signal émis
sigma_ne_K_codes_RS = sqrt((R_K_RS*P_signal_e_K_codes_RS * Ns)./(2*log2(M)*Eb_N02));
dim1 = length(Eb_N02);
dim2_RS = length(signal_e_K_codes_RS);
n_I_K_codes_RS = zeros(dim1, dim2_RS);
n_Q_K_codes_RS = zeros(dim1, dim2_RS);
n_e_K_codes_RS = zeros(dim1, dim2_RS);
for i=1:dim1
    n_I_K_codes_RS(i,:) = sigma_ne_K_codes_RS(i).*randn(1,dim2_RS);
    n_Q_K_codes_RS(i,:) = sigma_ne_K_codes_RS(i).*randn(1,dim2_RS);
    n_e_K_codes_RS(i,:) = n_I_K_codes_RS(i,:) + 1i* n_Q_K_codes_RS(i,:);
end

r_K_codes_RS = repmat(signal_e_K_codes_RS, dim1, 1) + n_e_K_codes_RS;

% Réception
z_K_codes = zeros(dim1, dim2);
z_K_codes_RS = zeros(dim1, dim2_RS);
for i=1:dim1
    z_K_codes(i,:) = filter(h_reception,1, r_K_codes(i,:));
    z_K_codes_RS(i,:) = filter(h_reception,1, r_K_codes_RS(i,:));
end
clear i;

% Echantillonnage
z_K_codes = z_K_codes(:, n0: Ns :end);
z_K_codes_RS = z_K_codes_RS(:, n0: Ns :end);

% Démodulation
d_K_codes = sign(real(z_K_codes)) + 1i * sign(imag(z_K_codes));
d_K_codes_RS = sign(real(z_K_codes_RS)) + 1i * sign(imag(z_K_codes_RS));

% Décision
bits_recus_K_codes = zeros(dim1, N_K);
bits_recus_K_codes_RS = zeros(dim1, N_K_RS);
for i = 1:dim1
    bits_recus_K_codes(i,:) = (qamdemod(d_K_codes(i,:).', M, 'OutputType', 'bit')).';
    bits_recus_K_codes_RS(i,:) = (qamdemod(d_K_codes_RS(i,:).', M, 'OutputType', 'bit')).';
end
clear ligne_ii;

% Décodage
bits_K_decodes = zeros(1, length(bits_K)); 
bits_K_decodes_RS = zeros(1, length(bits_K)); 
for i=1:dim1
    bits_K_decodes(i,:) = vitdec(bits_recus_K_codes(i,:),trellis,tblen,'trunc','hard', matrice_p);
    v = vitdec(bits_recus_K_codes_RS(i,:),trellis,tblen,'trunc','hard', matrice_p);
    v = deintrlv(v, vect_perm);
    v = step(decoder, v.');
    bits_K_decodes_RS(i,:) = v.'; 
end

% Calcul du taux d'erreur binaire
    % Sans entrelaceur
TEB_K = zeros(1,dim1);
for i=1:dim1
    TEB_K(i) = mean(abs(bits_K-bits_K_decodes(i,:)));
end
clear i;
    % Avec entrelaceur
TEB_K_RS = zeros(1,dim1);
for i=1:dim1
     TEB_K_RS(i) = mean(abs(bits_K-bits_K_decodes_RS(i,:)));
end
clear i;

% Comparaison entre TEB avec entrelaceur et TEB sans entrelaceur
figure()
semilogy(Eb_N0_dB2, TEB_K);
hold on
grid on
semilogy(Eb_N0_dB2, TEB_K_RS);
legend("Sans entrelaceur","Avec entrelaceur")
xlabel('E_b/N_0 (dB)')
ylabel('TEB')
title('Comparaison entre TEB avec entrelaceur et TEB sans entrelaceur ')