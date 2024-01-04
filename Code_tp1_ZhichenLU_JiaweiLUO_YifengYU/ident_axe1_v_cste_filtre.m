%% CE SCRIPT PERMET D'IDENTIFIER LES PARAMETRES DU MODELE
%% DE L'AXE 2 POUR DES MOUVEMENTS A VITESSE CONSTANTE
%% v1 G. MOREL - 2005/12/29
%% v2 V. PADOIS - 2007/10/31
%% v3 V. PADOIS - 2012/12/03

clear all; %% efface toutes les variables existantes
load releve_vit_cste_axe1; %% charge les releves experimentaux

%% Parametres connus a priori:
kc1=0.0525; %% constante de couple de l'axe 1.
N1=20.25; %% inverse du rapport de reduction de l'axe 1.

%% Construction de la matrice Y.
for(i=1:length(ifil1))  %% length(i1) = nombre d'echantillons.
    Y(i,1:4) = [ cos(q1(i)) sign(qpfil1(i)) qpfil1(i) 1 ];
    u(i,1) = kc1*N1*ifil1(i);
end
%% Calcul des parametres
p1=pinv(Y)*u;

%% Affichage des resultats.
format long
disp('Parametres estimes a partir des donnees brutes :');
p1'

figure(1)
clf; %% clear figure
h=plot3(q1,qpfil1,kc1*N1*ifil1,'m');
set(h,'LineWidth',0.5);
hold on; %% permet de conserver le graphique et d'en ajouter d'autres sur la meme fig.
h=plot3(q1,qpfil1,Y*p1,'r');
set(h,'LineWidth',1.5);
title('Resultats de l''identification avec filtrage');
legend('\Gamma_1 filtre', 'modele');
grid on;
