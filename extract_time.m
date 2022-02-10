

function [sati,minuti]=extract_time(I)
%Ime funkcije:extract_time
%Funkcija se koristi za racunanje vremena,odnosno prepoznavanje trenutnog
%vremena na osnovu slika zadatih satova.Zadatak ove funkcije je da prepozna
%kazaljke sata,odnosno kazaljku za sate i minute,ako one obe postoje i da
%na osnovu prepoznatih kazaljki vrati odgovarajuce vreme.Funkcija radi i za
%slucajeve kada imamo manje od 2 kazaljke,odnosno za slucaj kada su kazaljke
%preklopljene,pri cemu vraca odgovarajuce resenje. Funkciji se prosledjuje
%slika u odgovarajucem formatu i na njoj se vrsi dalja obrada.Koristimo
%Hafovu transformaciju kako bismo izvrsili detekciju ivica i dosli do
%odgovarajucih kazaljki,a uz pomoc odredjenih celija date strukture kao sto
%su ugao i duzina,dolazimo do saznanja koja je koja kazaljka i pod kojim se
%uglom nalazi tako da to mozemo konvertovati u odgovarajuci broj,odnosno
%sate i minute.
%
%Izgled funkcije:
%
%[sati,minuti]=extract_time(I); Ulazna slika I je u odgovarajucem formatu,a
%izlaz funkcije su odgovarajuci brojevi koji oznacavaju sate i minute.
%
%Funkcija nema podrazumevane vrednosti.Funkciji se eksplicitno mora
%proslediti slika u datom formatu.
%
%Primer:
%
%     I1=im2double(rgb2gray(imread('clock1.png')));
%     [sati1,minuti1]=extract_time(I1);
%     disp(['Prvi sat: ', num2str(sati1), ' :  ' ,num2str(minuti1)]);
%
%
%
%See also: extract_time_bonus
%
%Funkcija extract_time_bonus predstavlja nadogradnju date funkcije koja
%vraca i informaciju u trenutnom broju sekundi.
%
%
% Dan kreacije: 28.12.2019. (Petkovic Uros)
% Poslednje izmene: 28.12.2019. (Petkovic Uros)
%

I=medfilt2(I,[7 7]);     %Radi bolje obrade sliku prvo prepustamo kroz median filtar
[M,N]=size(I);           %Nalazimo dimenzije matrice
Z=M/9; C=M/7;            %Karakteristicne tacke za uporedjivanje praga piksela
[E,~]=edge(I,'canny',0.2,2);    %Nalazimo ivice
[H, T, R] = hough(E);           %Koristimo Hafovu transformaciju
P = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));  %Nalazimo pikove
lines = houghlines(E, T, R, P,'FillGap',25,'MinLength',80); %Nalazimo date linije
centar=[round(N/2) round(M/2)];          %Nalazimo centar slike koji ce biti znacajan za dalju obradu
%Centar je N/2 i M/2 iz razloga sto je obrnuto kad je u pitanju slika i
%obicna matrica
brojac=[];          %Postavljamo brojac na 0
flag=0;             %Postavljamo flag na 0

%Posto nam funkcije houghlines vraca vise od 3 linije,u sledecem foru
%prolazimo kroz sve linije,pri cemu cemo izabrati samo one 3 linije koje
%predstavljaju nase kazaljke,ako imamo dve linije koje su suvise
%blizu,uzecemo samo jednu,npr. slucaj sa satom 9. koji vraca dve linije za
%istu kazaljku sa obe strane,takodje,imamo i dodatne linije koje uopste ne
%predstavljaju kazaljke,pa ih treba odbaciti.Niz brojac nam vraca one
%vrednosti na kojima se nalaze nase kazaljke koje treba uzeti iz strukture
%lines
for k=1:length(lines)   %Prolazimo kroz ceo niz linija
    flag=0;
    xy=[lines(k).point1;lines(k).point2];    %Izdvajamo odgovarajuce pocetnu i krajnju koordinatu linije 
    if(sqrt((xy(1,1)-centar(1))^2+(xy(1,2)-centar(2))^2)<Z ||...
            sqrt((xy(2,1)-centar(1))^2+(xy(2,2)-centar(2))^2)<Z) %Uporedjujemo sa datim pragom blizine
        if(sqrt((xy(1,1)-xy(2,1))^2+(xy(1,2)-xy(2,2))^2)>C)
            if (~isempty(brojac))
                for i=1:length(brojac)
                    if(abs(lines(k).theta-lines(brojac(i)).theta)<5)
                        flag=1;       %Ako smo je nasli,postavljamo flag na 1
                    end
                end
                if(~flag)
                    brojac=[brojac k];  %Dodajemo na vec postojeci niz brojaca odg. vrednost indeksa linije
                end
            else
                brojac=[brojac k];
            end
        end
    end
    if (k>5)       %Ne treba nam vise od 5 linija,stoga izlazimo iz fora
        break
    end
end
linije=lines;                    %Pravim novi niz linije koji ce biti duzine 1,2 ili 3 u zavisnosti
if (length(brojac)==1)           %od toga koji broj kazaljki imamo
    linije(1)=lines(brojac);     %Ako imamo jednu kazaljku,smestamo je u linije
end
if (length(brojac)==2)           %Ako imamo dve,smestamo ih u prva dva mesta
    linije(1)=lines(brojac(1));
    linije(2)=lines(brojac(2));
end
if (length(brojac)==3)           %Ako imamo tri,smestamo ih u prva tri mesta
    linije(1)=lines(brojac(1));
    linije(2)=lines(brojac(2));  
    linije(3)=lines(brojac(3));
end
linije=linije(1:length(brojac)); %Na samom kraju niz linije ogranicavamo na onu duzinu koliko
%u stvari imamo kazaljki,to ce biti 1,2 ili 3,odbacujemo ostatak linija
%koje su nam visak i nisu nam potrebne za dalju obradu

Hx = [-1 -2 -1; 0 0 0; 1 2 1];    %Matrice Sobelovog operatora
Hy = Hx';                %Vertikalni i horiznotalni gradijenti
Gx = imfilter(I, Hx, 'replicate', 'same');
Gy = imfilter(I, Hy, 'replicate', 'same');
Gm = (Gx.^2 + Gy.^2).^0.5;   %Magnituda gradijenta
Ga = atan(Gy./Gx);           %Ugao gradijenta
%Racunamo magnitudu gradijenta zato sto ce nam biti potrebna kasnije u
%nadogradnji funkcije kada cemo uporedjivati sekundaru sa ostale dve
%kazaljke

for i=1:length(linije)               %Racunam duzinu za svaku kazaljku i magnitudu
    xy=[linije(i).point1;linije(i).point2];
    linije(i).duzina=sqrt((xy(1,1)-xy(2,1))^2+(xy(1,2)-xy(2,2))^2);
    magx=round((xy(1,1)+xy(2,1))/2); %Koristim tacku na sredini linije zbog krajeva koji su losiji
    magy=round((xy(1,2)+xy(2,2))/2); %Kazaljka moze dodirivati broj i onda samim tim imati vecu
    linije(i).mag=Gm(magy,magx);     %magnitudu koja ce dovoditi do losijih rezultata
end        
%Koristeci informaciju o magnitudi i smatrajuci da sekundara ima najmanju
%magnitudu,poredimo sve tri linije i smestamo ih tako da je treca ona sa
%najmanjom magnitudom,odnosno sekundara tipicnim algoritmom za sortiranje
if (length(linije)==3)     
       if (linije(1).mag<linije(2).mag)
       pom=linije(1);
       linije(1)=linije(2);
       linije(2)=pom;
       end
       if (linije(2).mag<linije(3).mag)
       pom=linije(2);
       linije(2)=linije(3);
       linije(3)=pom;
       end
end
%Sledeca ideja je da poredjam sve kazaljke tako da je uvek prva ona za sate,druga za minute
%a treca za sekunde iz razloga jer ce mi ta informacija u daljem toku
%programa biti jako bitna,pa da stalno ne bih ispitivao u kodu i pravio
%dodatnu slozenost,odmah cu ih tako razvrstati i onda samo koristiti vec
%dobijenu informaciju za ilustraciju ostatka koda
if (length(linije)>=2)  %Namestam da je prva uvek satara,a druga minutara
   if (linije(1).duzina>linije(2).duzina)
       pom=linije(1);
       linije(1)=linije(2);
       linije(2)=pom;
   end
end  
%Takodje,posto cu gledati u kom se "kvadrantu" sata nalazim,bitno mi je da
%znam koja mi je koja tacka,a obzirom da jedna tacka mora biti u centru
%sata,postavljam uvek da mi je to prva tacka,a druga tacka predstavlja onu
%sa drugog kraja koju cu koristiti u daljoj obradi,takodje.

%Namestam da je uvek prva tacka ona bliza centru za sataru
xy=[linije(1).point1;linije(1).point2];
if ((xy(2,1)-centar(1))<50) && ((xy(2,2)-centar(2))<80)
   pom=linije(1).point1;
   linije(1).point1=linije(1).point2;
   linije(1).point2=pom;
end
%Isto radim za minutaru ako je imam
if (length(linije)>=2)
   xy=[linije(2).point1;linije(2).point2];
   if ((xy(2,1)-centar(1))<50) && ((xy(2,2)-centar(2))<80)
      pom=linije(2).point1;
      linije(2).point1=linije(2).point2;
      linije(2).point2=pom;
   end    
end
%Isto to za sekundaru ako je imam
if (length(linije)==3)
   xy=[linije(3).point1;linije(3).point2];
   if ((xy(2,1)-centar(1))<50) && ((xy(2,2)-centar(2))<80)
      pom=linije(3).point1;
      linije(3).point1=linije(3).point2;
      linije(3).point2=pom;
   end    
end

%NAKON OVOGA SU POREDJANE SVE KAZALJKE REDOM SATARA,MINUTARA,SEKUNDARA I TO
%TAKO DA JE PRVA TACKA UVEK ONA BLIZA CENTRU

%Racunanje sati i minuta na osnovu kazaljki

if (length(linije)==1)          %Ako imam samo jednu kazaljku,odnosno satara i minutara su se preklopile
    if (linije(1).theta<0)      %Ako je ugao negativan,oduzimam njegobu apsolutnu vrednost od 360 stepeni
       UgaoSata=360-abs(linije(1).theta);
       UgaoMin=360-abs(linije(1).theta);
    else
       UgaoSata=linije(1).theta;   %Ako je ugao pozitivan,ostaje onakav kakav je bio
       UgaoMin=linije(1).theta;
    end 
else
   if (linije(1).theta<0)          %Za slucaj sa dve kazaljke radimo isto i za sate i za minute
       UgaoSata=360-abs(linije(1).theta);
   else
       UgaoSata=linije(1).theta;
   end  

   if (linije(1).point2(2)>centar(2))  %Ako je y komponenta druge kranje tacke veca od y komponente centra
       UgaoSata=mod((UgaoSata+180), 360);%znaci da se nalazim u donjoj polovini slike,u suprotnom u gornjoj
   end                                   %a ako sam u donjoj,onda moram dodati 180 stepeni,pri cemu radim mod
                                         %ako slucajno predjem 360 stepeni
                                         %da nemam vise od njega
   if (linije(2).theta<0) 
       UgaoMin=360-abs(linije(2).theta);
   else
       UgaoMin=linije(2).theta;
   end 

   if (linije(2).point2(2)>centar(2))    %Isto vazi i za minute
       UgaoMin=mod((UgaoMin+180), 360);
   end
end

sati=fix(UgaoSata/30);  %Posto je moguce imati 12 sati na raspolaganju,ugao delim sa 30 i koristim
if (sati==0)            %Funkciju fix koja zaokruzuje na manji ceo broj jer nam tako odgovara za sate
    sati=12;            %Ako je nula sati,napisi da je 12 sati
end
minuti=round(UgaoMin/6);%Takodje,mogucih vrednosti za minute je 60,stoga delimo sa 6,a mozemo koristiti
if (minuti==60)         %round jer je ovo ipak neka grublja procena minuta,pa ako je presla pola,moze se
    minuti=0;           %reci kao da je to skoro vec sledeci minut,sto za sate ne bi bio slucaj,jer bi
end                     %tada greska bila jako velika
                        %Ako je 60 minuta,napisi da je nula minuta
  
%U sledecoj funkciji,bonus funkciji,ostavljen je deo koda koji za sve satove pokrece slike koje ce                      
%biti prikazane u izvestaju domaceg zadatka                        
end
