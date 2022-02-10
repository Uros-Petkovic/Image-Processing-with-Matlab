
%% Prva tacka a)nacin

clc; clear all;

f = im2double(imread('girl_ht.tif'));     %Citamo sliku i vrsimo njen prikaz
figure; imshow(f);
[M, N] = size(f);
P = 2*M-1; Q = 2*N-1;                     %Prosirujemo granice zbog Furijeove transformacije
Fp = fftshift(fft2(f, P, Q));             %Vrsimo Furijeovu transformaciju nad slikom
H = lpfilter1('gaussian', P, Q, 65);       %Pravimo low pass filtar
H = fftshift(H);                          %Prikaz spektra filtra
figure; imshow(log(1+abs(H)), []);
figure; imshow(log(1+abs(Fp)), []);       %Prikaz spektra ulazne slike
Gp = Fp.*H;
figure; imshow(log(1+abs(Gp)), []);       %Prikaz spektra izlazne slike
gp = ifft2(ifftshift(Gp));    
g = gp(1:M, 1:N);
figure; imshow(g);                        %Prikaz izlazne slike nakon lp filtra

gn = log(1 + 1*g)/log(2);                           %Kompresija kontrasta log funkcijom
figure; imshow(gn);

gn=adapthisteq(gn,'ClipLimit', 0.0005,'NumTiles', [30 30]);  %Adaptivna ekvalizacija za dodatno poboljsanje
figure; imshow(gn);                                          %Prikaz slike
 
%Ova slika se moze dobiti pravljenjem elipsoidnog Gausovog lowpass filtra,
%pri cemu se odmah dobija dobra slika,samo malo tamna i mutnija,pa iz
%subjektivnih razloga radimo kompresiju kontrasta i dodatnu adaptivnu
%ekvalizaciju,a ista slika uradjena je i na drugi nacin koriscenjem kruznog
%Gausovog lowpass filtra pri cemu se javljaju dodatne horizontalne linije
%nakon obrade slike i potrebno je dodatno filtrirati sliku vertikalnim
%filtrom

%Propustanjem slike kroz lowpass filtar propustamo samo sredisnju komponentu spektra
%dok ostale komponente odbacujemo,ovaj postupak se mogao mozda i uraditi na
%neki drugi nacin,koriscenjem notch filtra sa matricom C koja sadrzi sve
%komponente koje zelimo odstraniti,ali uvideo sam da je matrica velika i
%samo izvrsavanje ovog dela zadatka je trajalo previse,a i nije davalo
%zadovoljavajuce rezultate kao ovaj metod,medjutim,ovde nije kraj
%obrade,nakon filtriranja se na slici pojavljuju horizontalne crne
%linije,pa dalje vrsimo obradu slike uz pomoc vertikalnog filtra

%% Prva tacka b)nacin

clc; clear all;

f = im2double(imread('girl_ht.tif'));     %Citamo sliku i vrsimo njen prikaz
figure; imshow(f);
[M, N] = size(f);
P = 2*M-1; Q = 2*N-1;                     %Prosirujemo granice zbog Furijeove transformacije
Fp = fftshift(fft2(f, P, Q));             %Vrsimo Furijeovu transformaciju nad slikom
H = lpfilter('gaussian', P, Q, 65);       %Pravimo low pass filtar
H = fftshift(H);                          %Prikaz spektra filtra
figure; imshow(log(1+abs(H)), []);
figure; imshow(log(1+abs(Fp)), []);       %Prikaz spektra ulazne slike
Gp = Fp.*H;
figure; imshow(log(1+abs(Gp)), []);       %Prikaz spektra izlazne slike
gp = ifft2(ifftshift(Gp));    
g = gp(1:M, 1:N);
figure; imshow(g);                        %Prikaz izlazne slike nakon lp filtra


[M, N] = size(g);                                    %Trazimo dimenzije slike
F = fft2(g, M, N);                                   %Vrsimo Furijeovu transformaciju slike
figure; imshow(log(1 + abs(fftshift(F))),[]);        %Prikazujemo spektar ulazne slike
H = recnotch('reject', 'vertical', M, N,11,70,70);   %Vrsimo uklanjanje horizontalnih linija
figure; imshow(log(1 + abs(fftshift(H))),[]);        %Prikaz spektra filtra
G = H.*F;                                            %Filtriranje slike
figure; imshow(log(1 + abs(fftshift(G))),[]);        %Prikaz spektra filtrirane slike
g1 = ifft2(G);                                       %Inverzija Furijeove transformacije
figure; imshow(g1);                                  %Prikaz izlazne slike

%Kod je pokretan i za manje parametre vertikalnog filtra gde je uoceno da vertikalne linije
%nisu u potpunosti uklonjene,pa sam ostavio date parametre za koje sam
%dobio dobre rezultate,a prikaz slike i za losiji rezultat dacu u izvestaju
%domaceg zadatka,sto se tice odabira parametara,na pamet mi je palo da
%stavim parametar 70 zato sto sam gledajuci spektar prvobitne slike uvideo
%da se elementi spektra javljaju na svakih 70 odbiraka,sto je u stvari i
%rezultovalo najboljem otklanjanju novonastalih horizontalnih linija

%Dalja obrada slike je cisto subjektivna,posto je nakon prvobitne obrade
%slike novonastala slika malo mutnija i tamnija,prvo cemo izvrsiti
%kompresiju kontrasta log funkcijom,a potom malo adaptivne ekvalizacije
%radi isticanja ivica
                                                     
gn = log(1 + 1*g1)/log(2);                                  %Kompresija kontrasta log funkcijom
figure; imshow(gn);

gn=adapthisteq(gn,'ClipLimit', 0.001,'NumTiles', [45 34]);  %Adaptivna ekvalizacija za dodatno poboljsanje
figure; imshow(gn);                                         %Prikaz slike

%% Druga tacka

clc; clear all;

Ibezsuma=im2double(imread('lena.tif'));
i=im2double(imread('lena_noise.tif'));    %Citanje i prikaz ulazne slike
noise=Ibezsuma-i;
varsuma=var(noise(:));      %Provereno da se treba dobiti vrednost 0.007
figure; imshow(i);

%Ideja je da za svaku dimenziju prozora odredimo matricu lokalnih varijansi
%i da nakon toga prikazemo histogram lokalnih varijansi,na tom histogramu
%vidimo koliko odbiraka matrice ima koju varijansu,a posto su pikseli
%priblizno uniformno raspodeljeni,mozemo naci glavnu varijansu tako sto
%cemo naci pik histograma i procitati vrednost na horizontalnoj osi koja u
%stvari predstavlja vrednosti datih varijansi,ideja je oduzeti matricu gvar^2
%od matrice (gvar)^2,ali postoji i prostiji i manje zametan nacin,postoji
%funkcija stdfilt koja za zadatu ulaznu sliku i sirinu prozora vraca
%matricu lokalnih standardnih devijacija slike,a posto nama trebaju lokalne
%varijanse,samo cemo tu matricu dici na kvadrat i dobiti zeljenu
%matricu,potom prikazati histogram i procitati vrednosti glavne varijanse

j=1;
for k=3:6:27
       ksus=ones(k,k);                  %Biranje prozora od 3x3 do 27x27
       I=stdfilt(i,ksus);               %Funkcija koja vraca matricu lokalnih standardnih devijacija
       J=I.^2;                          %Matrica lokalnih varijansi
       figure(j);histogram(J);          %Histogram varijansi
       xlabel('Varijansa suma');
       ylabel('Broj odbiraka sa datom varijansom');
       title(['Histogram lokalnih varijansi za k = ', num2str(k)]);
       j=j+1; 
    
end
% Za vece vrednosti susedstva dobija se tacnija vrednost varijanse suma ako
%gledamo pocetak diskretnog odbirka,a ne njegovu sredinu,nakon K=15,pocetak
%najveceg odbirka je uvek u 0.007,sto odgovara nasoj varijansi suma,ali 
%histogram najvise lici na Gausovu raspodelu za manja susedstva,sto je
%susedstvo vece,histogram je deformisaniji,nalazimo za svaki histogram
%komponentu na kojoj se nalazi maksimalni broj elemenata matrice J i ta
%komponenta predstavlja varijansu suma. Za date vrednosti susedstva
%dobijene su sledece vrednosti:

%k=3     var=0.0057
%k=9     var=0.0067
%k=15    var=0.007   priblizno za k=15 najbolji rezultat
%k=21    var=0.007
%k=27    var=0.007

%Funkcijom roipoly je, takodje, prethodno utvrdjeno da je varijansa jednaka
%var=0.007, kao sto smo to potvrdili i racunanjem varijanse suma na pocetku
%sto se poklapa sa tvrdjenjem da je neko optimalno k za procenu varijanse
%recimo za k=15 i za svako vece k
%% Treca tacka

%Ideja je otkloniti smetnje nastale na slici prilikom pokreta.Posto su nam
%dati podaci o tome kako se aparat kretao prilikom snimanja slike,mozemo
%iskoristiti te podatke kao i dobijenu sliku da izvrsimo otklanjanje
%smetnje tako sto cemo upotrebiti Vinerov filtar,a pre toga uraditi
%Furijeovu transformaciju slike i pokreta kamere. Posto pokret kamere nije
%istih dimenzija kao slika,u funkciji fft2 prosledjujemo i parametre
%dimenzija na kojima hocemo da uradimo Furijeovu transformaciju pokreta
%kamere,za razlicite parametre Vinerovog filtra dobijamo razlicite
%rezultate,kada je slobodan clan u imeniocu Vinerovog filtra manji,slika ce
%biti sve ostrija,ali doci ce do pojave belicastih tacaka na slici,to moze biti
%i uzrok zakucavanja visokih vrednosti na 1 i niskih na 0,dok ce
%za manje parametre,slika biti manje ostra i lepsa i jasnija,stoga nece
%biti potrebna neka preterana dodatna obrada slike.U izvestaju bice
%prikazana oba slucaja radi uvida u problem


clc; clear all;

gg=imread('etf_blur.tif');                     %Citamo sliku
k=im2double(imread('kernel.tif'));             %Citamo pokret kamere
e=im2double(gg);
figure; imshow(e);                             %Prikaz slike i pokreta kamere
figure; imshow(k);
[M, N] = size(e);
G = fftshift(fft2(e));                         %Furijeova transformacia slike
K = fftshift(fft2(k,M,N));                     %Prosirena F.t. pokreta kamere na dimenzije slike
W = (abs(K).^2)./(abs(K).^2 + 3);              %Vinerov filtar,koristimo veci broj u imeniocu
Fest = (G./K).*W;                              %Estimacija slike
figure; imshow(log(1+abs(G)), []);         
figure; imshow(log(1+abs(Fest)), []);          %Prikaz spektra estimirane slike
fest = ifft2(ifftshift(Fest));                 
fest(fest<0)=0;                                %Negativne vrednosti zakucam na 0
fest=fest./mean(fest(:)).*mean(e(:));          %Delim sa sr. vr. izlazne slike i mnozim sa sr. vr. ulazne slike
fest=fest(1:530,1:800);                        %Odsecanje crnog prozora,cisto iz subjektivnog razloga
figure; imshow(fest);                          %Prikaz estimirane slike


%% Cetvrta tacka

%Potrebno je izvrsiti nelokalno usrednjavanje ulazne zasumljene slike
%Kao varijansa suma slike uzima se rezultat iz druge tacke koji iznosi
%var=0.0075.Kod pokrecemo za razlicite vrednosti prozora susedstva i
%oblasti pretrazivanja za K=3,5,9 i S=15,33,51,dobijeni rezultati bice
%prikazani u izvestaju domaceg zadatka,kodovi su pokretani za vrednost
%faktora odsumljavanja 0.05 kao sa predavanja,iako smo mogli koristiti i
%druge parametre,sto je faktor odsumljavanja veci i stepen odsumljavanja ce
%biti veci,ali uzeli smo neku zlatnu sredinu,sto se tice nelokalnog
%usrednjavanja,najvise iteracija ce biti za kombinaciju 5x5 i 51x51 koja se
%i najvise izvrsavala,a optimalno resenje u skladu sa dobijenim rezultatima
%i vremenom izvrsavanja bice izreceno u izvestaju,zbog potrebe vise
%tacaka,na kraju je kod pokretan i za vrednosti K=7 i S=25

close all;
clc;
clear all;
 
Ibezsuma=im2double(imread('lena.tif'));
I = im2double(imread('lena_noise.tif'));    %Citanje ulazne slike
var=0.007;                                  %Procenjena varijansa iz drugog zadatka
h=0.05;                                     %Faktor odsumljivanja uzet proizvoljno

K=9;                              %Za svaku kombinaciju parametara pozivamo funkciju i prikazujemo
S=15;                             %dobijene rezultate na istom grafiku radi uporedjivanja
                                  %Ja sam ovde izabrao samo jednu kombinaciju,ali moze se pokrenuti za sve 
 
start=tic;                         %Oznacava pocetak odbrojavanja
J=dos_non_local_means(I,K,S,var,h);    %Izlazna slika iz funkcije

figure;                         %Prikazujemo svaku izlaznu sliku na zajednickom grafiku sa ulaznom slikom
subplot(1,2,1),imshow(I),title('Zasumljena slika'),colormap(gray);
subplot(1,2,2),imshow(J),title('Izlazna slika'),colormap(gray);
title(['Dobijeni rezultat za kombinaciju parametara K= ',num2str(K),' i S= ',num2str(S)]);
tElapsed = toc(start)/60;       %Racunamo potrebno vreme za izvrsenje programa za svaku sliku
noise=Ibezsuma-J;            %Mogli smo i var(noise(:)),ali imamo ugradjenu funkciju
PSNR=psnr(J,Ibezsuma);       %Racunamo odnos signal/sum za svaku sliku

%Mislio sam prvo uz pomoc suma,zato sam i racunao noise,ali nadjoh posle
%ugradjenu funkciju psnr,pa sam nju iskoristio


%Ovde smo pokretali za svaku kombinaciju parametara kod,a potom smo izabrali 
%najbolju kombinaciju cije sam parametre ostavio,a prikazane rezultate svih slika stavio u
%izvestaj domaceg zadatka,stoga se ovaj kod moze pokrenuti za bilo koju
%drugu kombinaciju parametara ako korisnik ima zelju za tim,ovo se moglo
%implementirati kroz dva fora koji prolaze kroz sve kombinacije,ali bi
%izvrsavanje tog koda bilo jako dugo,pa sam se zadrzao na ovom resenju.

                                                   %Petkovic Uros


%PSNR za 5 i 15   31.4060   vreme 4.4278 u minutima
%PSNR za 9 i 15   31.5950   vreme 4.6815 u minutima
%PSNR za 3 i 15   30.6792   vreme 4.3317 u minutima
%PSNR za 9 i 33   31.0513   vreme 22.1682 u minutima
%PSNR za 3 i 33   30.3309   vreme 19.4498 u minutima
%PSNR za 5 i 33   31.0961   vreme 19.9600 u minutima
%PSNR za 3 i 51   30.0215   vreme 45.7611 u minutima
%PSNR za 5 i 51   30.8041   vreme 46.4846 u minutima
%PSNR za 9 i 51   30.6520   vreme 51.0395 u minutima

%Dole pokrecem i za vrednosti sa 7 i 25


%% Pokretanje za sve slucajeve i plotovanje grafika PSNR i Vreme na dnu

%SKROLOVATI DOSTA! Moglo je sve u foru,ali ovo mi je bilo lakse zbog
%pokretanja

close all;
clc;
clear all;
 
Ibezsuma=im2double(imread('lena.tif'));
I = im2double(imread('lena_noise.tif'));   
var=0.007;                                 
h=0.05;                                    

start=tic;                   
J1=dos_non_local_means(I,3,33,var,h);
tElapsed1 = toc(start)/60;
PSNR1=psnr(J1,Ibezsuma);

figure; subplot(1,2,1),imshow(I),title('Zasumljena slika'),colormap(gray);
subplot(1,2,2),imshow(J1),title('Izlazna slika'),colormap(gray);
title(['Dobijeni rezultat za kombinaciju parametara K= ',num2str(3),' i S= ',num2str(33)]);

start=tic;                   
J2=dos_non_local_means(I,5,33,var,h);
tElapsed2 = toc(start)/60;
PSNR2=psnr(J2,Ibezsuma);

figure; subplot(1,2,1),imshow(I),title('Zasumljena slika'),colormap(gray);
subplot(1,2,2),imshow(J2),title('Izlazna slika'),colormap(gray);
title(['Dobijeni rezultat za kombinaciju parametara K= ',num2str(5),' i S= ',num2str(33)]);

start=tic;                   
J3=dos_non_local_means(I,3,51,var,h);
tElapsed3 = toc(start)/60;
PSNR3=psnr(J3,Ibezsuma);

figure; subplot(1,2,1),imshow(I),title('Zasumljena slika'),colormap(gray);
subplot(1,2,2),imshow(J3),title('Izlazna slika'),colormap(gray);
title(['Dobijeni rezultat za kombinaciju parametara K= ',num2str(3),' i S= ',num2str(51)]);

start=tic;                   
J4=dos_non_local_means(I,5,51,var,h);
tElapsed4 = toc(start)/60;
PSNR4=psnr(J4,Ibezsuma);

figure; subplot(1,2,1),imshow(I),title('Zasumljena slika'),colormap(gray);
subplot(1,2,2),imshow(J4),title('Izlazna slika'),colormap(gray);
title(['Dobijeni rezultat za kombinaciju parametara K= ',num2str(5),' i S= ',num2str(51)]);

start=tic;                   
J5=dos_non_local_means(I,9,51,var,h);
tElapsed5 = toc(start)/60;
PSNR5=psnr(J5,Ibezsuma);

figure; subplot(1,2,1),imshow(I),title('Zasumljena slika'),colormap(gray);
subplot(1,2,2),imshow(J5),title('Izlazna slika'),colormap(gray);
title(['Dobijeni rezultat za kombinaciju parametara K= ',num2str(9),' i S= ',num2str(51)]);

start=tic;                   
J6=dos_non_local_means(I,9,33,var,h);
tElapsed6 = toc(start)/60;
PSNR6=psnr(J6,Ibezsuma);

figure; subplot(1,2,1),imshow(I),title('Zasumljena slika'),colormap(gray);
subplot(1,2,2),imshow(J6),title('Izlazna slika'),colormap(gray);
title(['Dobijeni rezultat za kombinaciju parametara K= ',num2str(9),' i S= ',num2str(33)]);

start=tic;                   
J7=dos_non_local_means(I,3,15,var,h);
tElapsed7 = toc(start)/60;
PSNR7=psnr(J7,Ibezsuma);

figure; subplot(1,2,1),imshow(I),title('Zasumljena slika'),colormap(gray);
subplot(1,2,2),imshow(J7),title('Izlazna slika'),colormap(gray);
title(['Dobijeni rezultat za kombinaciju parametara K= ',num2str(3),' i S= ',num2str(15)]);

start=tic;                   
J8=dos_non_local_means(I,5,15,var,h);
tElapsed8 = toc(start)/60;
PSNR8=psnr(J8,Ibezsuma);

figure; subplot(1,2,1),imshow(I),title('Zasumljena slika'),colormap(gray);
subplot(1,2,2),imshow(J8),title('Izlazna slika'),colormap(gray);
title(['Dobijeni rezultat za kombinaciju parametara K= ',num2str(5),' i S= ',num2str(15)]);

start=tic;                   
J9=dos_non_local_means(I,9,15,var,h);
tElapsed9 = toc(start)/60;
PSNR9=psnr(J9,Ibezsuma);

figure; subplot(1,2,1),imshow(I),title('Zasumljena slika'),colormap(gray);
subplot(1,2,2),imshow(J9),title('Izlazna slika'),colormap(gray);
title(['Dobijeni rezultat za kombinaciju parametara K= ',num2str(9),' i S= ',num2str(15)]);

start=tic;                   
J10=dos_non_local_means(I,3,25,var,h);
tElapsed10 = toc(start)/60;
PSNR10=psnr(J10,Ibezsuma);

figure; subplot(1,2,1),imshow(I),title('Zasumljena slika'),colormap(gray);
subplot(1,2,2),imshow(J10),title('Izlazna slika'),colormap(gray);
title(['Dobijeni rezultat za kombinaciju parametara K= ',num2str(3),' i S= ',num2str(25)]);

start=tic;                   
J11=dos_non_local_means(I,5,25,var,h);
tElapsed11 = toc(start)/60;
PSNR11=psnr(J11,Ibezsuma);

figure; subplot(1,2,1),imshow(I),title('Zasumljena slika'),colormap(gray);
subplot(1,2,2),imshow(J11),title('Izlazna slika'),colormap(gray);
title(['Dobijeni rezultat za kombinaciju parametara K= ',num2str(5),' i S= ',num2str(25)]);

start=tic;                   
J12=dos_non_local_means(I,7,15,var,h);
tElapsed12 = toc(start)/60;
PSNR12=psnr(J12,Ibezsuma);

figure; subplot(1,2,1),imshow(I),title('Zasumljena slika'),colormap(gray);
subplot(1,2,2),imshow(J12),title('Izlazna slika'),colormap(gray);
title(['Dobijeni rezultat za kombinaciju parametara K= ',num2str(7),' i S= ',num2str(15)]);

start=tic;                   
J13=dos_non_local_means(I,7,25,var,h);
tElapsed13 = toc(start)/60;
PSNR13=psnr(J13,Ibezsuma);

figure; subplot(1,2,1),imshow(I),title('Zasumljena slika'),colormap(gray);
subplot(1,2,2),imshow(J13),title('Izlazna slika'),colormap(gray);
title(['Dobijeni rezultat za kombinaciju parametara K= ',num2str(7),' i S= ',num2str(25)]);

start=tic;                   
J14=dos_non_local_means(I,7,33,var,h);
tElapsed14 = toc(start)/60;
PSNR14=psnr(J14,Ibezsuma);

figure; subplot(1,2,1),imshow(I),title('Zasumljena slika'),colormap(gray);
subplot(1,2,2),imshow(J14),title('Izlazna slika'),colormap(gray);
title(['Dobijeni rezultat za kombinaciju parametara K= ',num2str(7),' i S= ',num2str(33)]);

start=tic;                   
J15=dos_non_local_means(I,7,51,var,h);
tElapsed15 = toc(start)/60;
PSNR15=psnr(J15,Ibezsuma);

figure; subplot(1,2,1),imshow(I),title('Zasumljena slika'),colormap(gray);
subplot(1,2,2),imshow(J15),title('Izlazna slika'),colormap(gray);
title(['Dobijeni rezultat za kombinaciju parametara K= ',num2str(7),' i S= ',num2str(51)]);

start=tic;                   
J16=dos_non_local_means(I,9,25,var,h);
tElapsed16 = toc(start)/60;
PSNR16=psnr(J16,Ibezsuma);

figure; subplot(1,2,1),imshow(I),title('Zasumljena slika'),colormap(gray);
subplot(1,2,2),imshow(J16),title('Izlazna slika'),colormap(gray);
title(['Dobijeni rezultat za kombinaciju parametara K= ',num2str(9),' i S= ',num2str(25)]);

PSNR=[PSNR1 PSNR2 PSNR3 PSNR4 PSNR5 PSNR6 PSNR7 PSNR8 PSNR9 PSNR10...
    PSNR11 PSNR12 PSNR13 PSNR14 PSNR15 PSNR16];
Vreme=[tElapsed1 tElapsed2 tElapsed3 tElapsed4 tElapsed5 tElapsed6...
    tElapsed7 tElapsed8 tElapsed9 tElapsed10 tElapsed11 tElapsed12...
    tElapsed13 tElapsed14 tElapsed15 tElapsed16];

t=1:16;
figure; stem(t,PSNR,'r*');
title('PSNR za razlicite vrednosti parametara K i S');
ylabel('Odnos u decibelima');
xlabel('Broj iteracije');
figure; stem(t,Vreme,'b*');
title('Vreme izvrsavanja u minutima za razlicite vrednosti parametara K i S');
ylabel('Vreme u minutima');
xlabel('Broj iteracije');

%Najveci PSNR za 9 i 15 ako gledamo samo dati skup,ako gledamo moj skup,
%najveci PSNR je za 7 i 15

