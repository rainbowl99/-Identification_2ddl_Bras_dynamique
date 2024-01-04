%% CE SCRIPT PERMET D'IDENTIFIER LES PARAMETRES DU MODELE
%% DE L'AXE 2 POUR DES MOUVEMENTS A VITESSE CONSTANTE
%% v1 G. MOREL - 2005/12/29
%% v2 V. PADOIS - 2007/10/31
%% v3 V. PADOIS - 2012/12/03

clear all; %% efface toutes les variables existantes
load releve_vit_cste_axe2; %% charge les releves experimentaux

%% Parametres connus a priori:
kc2=0.0525; %% constante de couple de l'axe 2.
N2=4.5; %% inverse du rapport de reduction de l'axe 2.

%% Construction de la matrice Y.
for(i=1:length(ifil2)) , %% length(i2) = nombre d'echantillons.
    Y(i,1:4) = [ cos(q2(i)) sign(qpfil2(i)) qpfil2(i) 1 ];
    u(i,1) = kc2*N2*ifil2(i);
end
%% Calcul des parametres
p2=pinv(Y)*u;

%% Affichage des resultats.
format long
disp('Parametres estimes a partir des donnees brutes :');
p2'

figure(2)
clf; %% clear figure
h=plot3(q2,qpfil2,kc2*N2*ifil2,'m');
set(h,'LineWidth',0.5);
hold on; %% permet de conserver le graphique et d'en ajouter d'autres sur la meme fig.
h=plot3(q2,qpfil2,Y*p2,'r');
set(h,'LineWidth',1.5);
title('Resultats de l''identification avec filtrage');
legend('\Gamma_2 filtre', 'modele');
grid on;
