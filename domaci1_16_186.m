
%% Zadatak 1 - Popravka kvaliteta slike

%Date slike se nalaze u istom folderu kao i .m fajl,stoga slike mogu
%procitati samim njihovim nazivom,posto se slike ne salju pri predaji
%domaceg,sekvenca dolaska do slike koja se treba procitati se naknadno moze
%promeniti na samoj odbrani domaceg zadatka. Petkovic Uros

clc;                                %Brisanje svega na pocetku
clear all;

Ibr=hdrread('bristolb.hdr');        %Citanje slike 'bristolb'
figure;
imshow(Ibr);                        %Prikaz ulazne slike 'bristolb'
Idis=imread('disney.png');          %Citanje slike 'disney'
figure;
imshow(Idis);                       %Prikaz ulazne slike 'disney'
Igir=imread('giraff.jpg');          %Citanje slike 'giraff'
figure;
imshow(Igir);                       %Prikaz ulazne slike 'giraff'
Ien=imread('enigma.png');           %Citanje slike 'enigma'
figure;
imshow(Ien);                        %Prikaz ulazne slike 'enigma'

s=12.92;
gama=1/2.4;                         %Parametri za prebacivanje iz linearnih
t=0.003138;                         %u sRGB format uz pomof transfer funkcije
f=0.055;

IbrR=Ibr(:,:,1);                    %Razdvojimo svaku komponentu posebno radi
IbrG=Ibr(:,:,2);                    %brze obrade
IbrB=Ibr(:,:,3);

N=numel(IbrR);                      %Ukupan broj elemenata u svakoj komponenti

for i=1:N                           %U zavisnosti od vrednosti ulazne komponente,
    if (IbrR(i)>t)                  %koristimo odgovarajucu formulu
       IbrR(i)=IbrR(i)^gama*(1+f)-f;
    else                     
        IbrR(i)=IbrR(i)*s;          %Dobijamo odgovarajuce R komponente
    end
end

for i=1:N
    if (IbrG(i)>t)
       IbrG(i)=IbrG(i)^gama*(1+f)-f;
    else
        IbrG(i)=IbrG(i)*s;          %Dobijamo odgovarajuce G komponente
    end
end

for i=1:N
    if (IbrB(i)>t)
       IbrB(i)=IbrB(i)^gama*(1+f)-f;
    else
        IbrB(i)=IbrB(i)*s;          %Dobijamo odgovarajuce B komponente
    end
end

Ibr(:,:,1)=IbrR;                    %Sve to pakujemo u novodobijenu sliku
Ibr(:,:,2)=IbrG;
Ibr(:,:,3)=IbrB;
Ibr=double(Ibr);                    %Prebacujemo u double format
figure;
imshow(Ibr);                        %Prikazujemo dobijenu sliku

Ibrgray=im2double(rgb2gray(Ibr));              %Siva 'bristolb'
Idisgray=im2double(rgb2gray(Idis));            %Siva 'disney'
figure; imshow(Ibrgray);                       %Prikaz slika
figure; imshow(Idisgray);

h=imhist(Ibrgray);                             %Crtamo histogram ulaznih slika                  
hnorm = h./numel(Ibrgray);                        
figure; bar(hnorm);                            
title('Normalizovani histogram ulaznog bristola');
h=imhist(Idisgray);                               
hnorm = h./numel(Idisgray);                        
figure; bar(hnorm);                           
title('Normalizovani histogram ulaznog disneya');

Jbrgray = Ibrgray.^0.38;                       %Radimo kompresiju sive slike 'bristolb'
figure; imshow(Jbrgray);                       %uz pomoc eksponencijalne funkcije
set(gcf, 'Name', 'Izlaz posle log funkcije za k = 1000');

w_sharp = [-1 -1 -1; -1 17 -1; -1 -1 -1]./9;   %Dodatno izostravamo sliku radi lepseg izgleda
Jbrgray =imfilter(Jbrgray, w_sharp, 'replicate');
figure; imshow(Jbrgray);

Jbrgray=adapthisteq(Jbrgray, 'ClipLimit', 0.007,'NumTiles', [45 34]);
figure; imshow(Jbrgray);                       %Radimo adaptivnu ekvalizaciju histograma za dodatno
                                               %podesavanje histograma
h=imhist(Jbrgray);                               
hnorm = h./numel(Jbrgray);                     %Izlazni histogram             
figure; bar(hnorm);                            
title('Normalizovani histogram izlaznog bristola');                                               
                                               
                                               
Jdisgray = log(1 + 100*Idisgray)/log(101);     %Radimo kompresiju kontrasta 'disney' slike koriscenjem
figure; imshow(Jdisgray);                      %log funkcije
set(gcf, 'Name', 'Izlaz posle log funkcije za k = 100');

Jdisgray=adapthisteq(Jdisgray, 'ClipLimit', 0.01,'NumTiles', [45 34]);
figure; imshow(Jdisgray);                      %Radimo adaptivnu ekv. histograma za dodatno ulepsavanje
                                               %Adaptivne ekv. histograma
                                               %nisu potrebne,to je izbor
                                               %urednika,sve je na
                                               %subjektivnom nivou
h=imhist(Jdisgray);                               
hnorm = h./numel(Jdisgray);                    %Izlazni histogram  
figure; bar(hnorm);                            
title('Normalizovani histogram izlaznog disneya'); 
                                                                                             
h=imhist(Igir);                                %Na osnovu histograma vidimo da su na slici zirafe
hnorm = h./numel(Igir);                        %komponente pretezno rasporedeljene od 90 do 160
figure; bar(hnorm);                            %vrednosti,pa cemo uzeti pocetnu i krajnju vrednost i
title('Normalizovani histogram slike');        %izvrsiti razvlacenje kontrasta tog opsega na opseg 0 do 255
Jgir = imadjust(Igir, [90 160]/255, [0 1],1.02);
figure; imshow(Jgir);
set(gcf, 'Name', 'Izlazna slika');

h=imhist(Jgir);                                %Izlazni histogram                     
hnorm = h./numel(Jgir);                        
figure; bar(hnorm);                            
title('Normalizovani histogram izlazne zirafe'); 

Ien=im2double(Ien);                            %Na datoj slici se vidi da imamo sum,zato koristimo median filtar
Jen = medfilt2(Ien, [11 11],'symmetric');      %koji sluzi za otklanjanje suma,medjutim,dobija se malo mutnija slika,
figure; imshow(Jen);                           %ali ulogu u tome igraju neki drugi artefakti slike,mozda smo i mogli
set(gcf, 'Name', 'Izlazna slika1');            %daljom obradom dobiti bolju sliku,ali zadrzali smo se na ovome


Idis=im2double(Idis);                          %Vrsimo konverziju u double format
figure; imshow(Idis);
set(gcf, 'Name', 'Ulazna slika');              %Prikazujemo ulaznu sliku

Ihsv = rgb2hsv(Idis);                          %Prebacujemo sliku u HSV kolor sistem
Ihsv(:,:,3) = Ihsv(:,:,3).^0.48;               %Podesavamo komponentu Value
Ihsv(:,:,2) = Ihsv(:,:,2).^0.92;               %Podesavamo komponentu Saturation
Jdis = hsv2rgb(Ihsv); figure; imshow(Jdis);    %Prikazujemo dobijenu sliku
set(gcf, 'Name', 'Slika obradjena u HSV');     %U ovom slucaju smo mogli vrsiti obradu i u YCbCr sistemu i u Lab sistemu
                                               %ali na uredniku je da
                                               %izabere ono sto mu se
                                               %najvise dopada

figure; imshow(Ibr);
set(gcf, 'Name', 'Ulazna slika');              %Prikaz ulazne slike

Ihsv = rgb2hsv(Ibr);                           %Takodje i u ovom delu radimo obradu u HSV sistemu iz licnog izbora
Ihsv(:,:,3) = Ihsv(:,:,3).^0.43;
Ihsv(:,:,2) = Ihsv(:,:,2).^0.83;
Jbr = hsv2rgb(Ihsv); figure; imshow(Jbr);
set(gcf, 'Name', 'Slika obradjena u HSV');     %Prikazujemo dobijenu sliku

%% Zadatak 2 - Testiranje funkcije dos_clhe za ekvalizaciju histograma

Iein=im2double(imread('einstein_lc.tif'));    %Ucitavamo sliku 'einstein_lc'

J1=dos_clhe(Iein,8,1);                     %Pozivamo implementiranu funkciju
figure; imshow(J1);                           %Prikazujemo dobijenu sliku

J2=dos_clhe(Ibr,8,1);                      %Pozivamo implementiranu funkciju
figure; imshow(J2);                           %Prikazujemo dobijenu sliku

