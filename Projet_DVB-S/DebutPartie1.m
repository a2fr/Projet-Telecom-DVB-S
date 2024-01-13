clear
close all

% Paramètres de l'étude
Rb = 3000;        % Débit binaire
Fe = 24000; % Fréquence d'échantillonnage
M = 4;         % Ordre de modulation
Rs = Rb / log2(M); % Débit symbole
Ns = sqrt(Fe / Rs);
alpha = 0.35;
L = 8;
Te = 1 / Fe;

% Transmission au format DVB-S2 avec mapping 8-PSK
% (Quadrature Phase Shift Keying) et filtre de mise en forme : racine de cosinus surélevé de roll off 0.2

% Génération du message binaire aléatoire
n = 20000;
message_binaire = randi([0, 1], 1, n);

% Mapping
symboles = qammod(message_binaire.', M, 'InputType', 'bit').';

% Suréchantillonnage
symboles_oversampled = kron(symboles, [1 zeros(1, Ns-1)]);
h = rcosdesign(alpha, L, Ns, 'sqrt');
signal_transmis = filter(h, 1, [symboles_oversampled zeros(1, length(h)-1)]);

% Canal AWGN passe-bas équivalent
Eb_N0_dB = 0:4;       % Valeurs pour le rapport signal sur bruit (Eb/No) en (dB) allant de 0 à 8 dB.
Eb_N0 = 10.^(Eb_N0_dB/10);     % Convertit les valeurs de Eb/No de dB en puissance réelle

puissance_signal = mean(abs(signal_transmis).^2);
ecart_type_bruit = sqrt((puissance_signal * Ns) ./ (2 * log2(M) * Eb_N0));

% Les dimensions pour les matrices de bruit.
nb_Eb_N0 = length(Eb_N0);
dim_signal_transmis = length(signal_transmis);

bruit_I = ecart_type_bruit' * randn(1, dim_signal_transmis);
bruit_Q = ecart_type_bruit' * randn(1, dim_signal_transmis);
bruit = bruit_I + 1i * bruit_Q;

signal_recu = repmat(signal_transmis, nb_Eb_N0, 1) + bruit;

% Réception
hr = h;
signal_filtre_recu = filter(hr, 1, signal_recu, [], 2);

% Réception sans bruit
signal_sans_bruit_recu = filter(h, 1, signal_transmis);

% Echantillonnage
n0 = length(h);
signal_recu_echantillonne = signal_filtre_recu(:, n0: Ns :end);
signal_sans_bruit_recu_echant = signal_sans_bruit_recu(n0: Ns :end);

% Décision
symboles_estimes = sign(real(signal_recu_echantillonne)) + 1i * sign(imag(signal_recu_echantillonne));
symboles_estimes_sans_bruit = sign(real(signal_sans_bruit_recu_echant)) + 1i * sign(imag(signal_sans_bruit_recu_echant));

bits_sans_bruit = qamdemod(symboles_estimes_sans_bruit.', M, 'OutputType', 'bit');
bits_sans_bruit = bits_sans_bruit';
bits_chapeau = zeros(nb_Eb_N0, length(message_binaire));

for i = 1:nb_Eb_N0
    ligne_i = qamdemod(symboles_estimes(i,:).', M, 'OutputType', 'bit');
    bits_chapeau(i,:) = ligne_i';
end

% Calcul du taux d'erreur binaire
TEB = mean(abs(message_binaire - bits_chapeau), 2).';

% TES sans bruit
TES_sans_bruit = mean(abs(symboles - symboles_estimes_sans_bruit));

% TES avec bruit
TES = mean(abs(symboles - symboles_estimes), 2).';

% Formules théoriques
TEB_theorique = qfunc(sqrt(2*Eb_N0));
TES_theorique = TEB_theorique*log2(M);

figure;
semilogy(Eb_N0_dB, TEB);
grid on;
hold on
semilogy(Eb_N0_dB, TEB_theorique);
grid on;
legend("TEB simulé",'TEB théorique')
xlabel('E_b/N_0 (dB)')
ylabel('TEB')
title('Comparaison entre TEB théorique et TEB simulé')

figure;
semilogy(Eb_N0_dB, TES);
grid on;
hold on
semilogy(Eb_N0_dB, TES_theorique);
grid on;
legend("TES simulé",'TES théorique')
xlabel('E_b/N_0 (dB)')
ylabel('TES')
title('Comparaison entre TES théorique et TES simulé')

% Tracé du diagramme de l'œil sans bruit
% eyediagram(signal_sans_bruit_recu, Ns*2);

% Tracé du diagramme de l'œil avec bruit
% eyediagram(signal_recu(nb_Eb_N0(end),:), Ns*2);

% Tracé de la réponse impulsionnelle
% figure;
% stem(h, 'b', 'LineWidth', 2);
% title('Réponse impulsionnelle du filtre en racine de cosinus surélevé');
% xlabel('Temps (échantillons)');
% ylabel('Amplitude');
% grid on;
