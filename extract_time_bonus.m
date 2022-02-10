

function [sati,minuti,sekunde]=extract_time_bonus(I)

%Ime funkcije:extract_time_bonus
%Funkcija se koristi za racunanje vremena,odnosno prepoznavanje trenutnog
%vremena na osnovu slika zadatih satova.Zadatak ove funkcije je da prepozna
%kazaljke sata,odnosno kazaljku za sate,minute i sekunde,ako one sve postoje i da
%na osnovu prepoznatih kazaljki vrati odgovarajuce vreme.Funkcija radi i za
%slucajeve kada imamo manje od 3 kazaljke,odnosno za slucaj kada su kazaljke
%preklopljene ili imamo samo kazaljku za sate i minute,pri cemu vraca
%odgovarajuce resenje. Funkciji se prosledjuje slika u odgovarajucem formatu
%i na njoj se vrsi dalja obrada.Koristimo Hafovu transformaciju kako bismo
%izvrsili detekciju ivica i dosli do odgovarajucih kazaljki,a uz pomoc 
%odredjenih celija date strukture kao sto su ugao i duzina,dolazimo do saznanja
%koja je koja kazaljka i pod kojim se uglom nalazi tako da to mozemo konvertovati
%u odgovarajuci broj,odnosno sate,minute i sekunde.
%
%Izgled funkcije:
%
%[sati,minuti,sekunde]=extract_time_bonus(I); Ulazna slika I je u odgovarajucem
%formatu,a izlaz funkcije su odgovarajuci brojevi koji oznacavaju
%sate,minute i sekunde.
%
%Ako se funkciji prosledi slika koja ima samo dve kazaljke ili cak jednu
%kazaljku,odnosno preklopljene su kazaljke sa sate i minute,funkcija kao
%informaciju o sekundama vraca -1.
%
%Primer:
%
%     I1=im2double(rgb2gray(imread('clock1.png')));
%     [sati1,minuti1,sekunde1]=extract_time_bonus(I1);
%     disp(['Prvi sat bonus: ', num2str(sati1), ' :  ' ,num2str(minuti1),' : ',num2str(sekunde1)]);
%
%
%
%See also: extract_time
%
%Funkcija extract_time_bonus predstavlja nadogradnju funkcije extract_time koja
%vraca informaciju u trenutnom broju sata i minuta.
%
%
% Dan kreacije: 28.12.2019. (Petkovic Uros)
% Poslednje izmene: 28.12.2019. (Petkovic Uros)
%




I1=I;              %Cuvamo jos jednom sliku pre obrade medijan filtrom jer cemo njom naci magnitudu
I=medfilt2(I,[5 5]);     %u nastavku programa,ovde,takodje,prolazimo kroz medijan filtar
[M,N]=size(I);           %Nalazimo dimenzije matrica
Z=M/9; C=M/7;            %Koristimo neke promenljive za pragove
[E,~]=edge(I,'canny');   %Nalazimo ivice
[H, T, R] = hough(E);    %Koristimo Hafovu transformaciju
P = houghpeaks(H,70,'threshold',ceil(0.1*max(H(:)))); %Nalazimo pikove
lines = houghlines(E, T, R, P,'FillGap',21,'MinLength',Z); %Nalazimo linije
centar=[round(N/2) round(M/2)];  %Nalazimo centar koji smo obrnuli iz istog razloga kao i u
%prethodnoj funkciji jer nije isto za matrice i za sliku

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
for k=1:length(lines)    %Prolazimo kroz ceo niz linija
    flag=0;
    xy=[lines(k).point1;lines(k).point2];   %Izdvajamo odgovarajuce pocetnu i krajnju koordinatu linije
    if(sqrt((xy(1,1)-centar(1))^2+(xy(1,2)-centar(2))^2)<Z ||...
            sqrt((xy(2,1)-centar(1))^2+(xy(2,2)-centar(2))^2)<Z)  %Uporedjujemo sa datim pragom blizine
        if(sqrt((xy(1,1)-xy(2,1))^2+(xy(1,2)-xy(2,2))^2)>C)
            if (~isempty(brojac))
                for i=1:length(brojac)
                    if(abs(lines(k).theta-lines(brojac(i)).theta)<5)
                        flag=1;     %Ako smo je nasli,postavljamo flag na 1
                    end
                end
                if(~flag)
                    brojac=[brojac k];
                end        %Dodajemo na vec postojeci niz brojaca odg. vrednost indeksa linije
            else
                brojac=[brojac k];
            end
        end
    end
    if (k>5)           %Ne treba nam vise od 5 linija,stoga izlazimo iz fora
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

if (length(linije)==3)     %Po magnitudi namestim da je treca uvek sekundara,najmanja magnituda
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

if (length(linije)==1)          %Ako imam jednu kazaljku,odnosno preklopile su se satara i minutara
    sekunde=-1;                 %ulazim u ovaj if,sekunde postavljam na -1 jer nemam sekundaru
    if (linije(1).theta<0)      %Ako je ugao manji od nula,prebacujem ga u odgovarajuci domen
       UgaoSata=360-abs(linije(1).theta);
       UgaoMin=360-abs(linije(1).theta);
    else 
       UgaoSata=linije(1).theta;%Ako je pozitivan,ostaje kakav je bio
       UgaoMin=linije(1).theta;
    end 
else
   sekunde=-1;                 %Slucaj sa dve kazaljke,isto postavljam sekunde na -1
   if (linije(1).theta<0)      %Negativan ugao prebacujem u opseg
       UgaoSata=360-abs(linije(1).theta);
   else
       UgaoSata=linije(1).theta;   %Pozitivan ostaje isti
   end  

   if (linije(1).point2(2)>centar(2))  %Nalazim se dole desno
       UgaoSata=mod((UgaoSata+180), 360);  %Dodajem 180 i radim mod da ne bih presao 360
   end

   if (linije(2).theta<0) 
       UgaoMin=360-abs(linije(2).theta);   %Negativan opet prebacujem u opseg
   else
       UgaoMin=linije(2).theta;            %Pozitivan ostaje isti
   end 

   if (linije(2).point2(2)>centar(2))      %Dole desno
       UgaoMin=mod((UgaoMin+180), 360);    %Prebacivanje u opseg
   end
end
if (length(linije)==3)                     %Slucaj sa 3 kazaljke
       if (linije(3).theta<0)              %Negativan ugao za sekunde prebacujem
          UgaoSekunde=360-abs(linije(3).theta);
       else
          UgaoSekunde=linije(3).theta;     %Ako je pozitivan,ostaje
       end  

       if (linije(3).point2(2)>centar(2))  %U kom se kvadrantu nalazim
          UgaoSekunde=mod((UgaoSekunde+180), 360);
       end
       sekunde=round(UgaoSekunde/6);       %Radim round jer mozemo zaokruziti i na vecu sekundu
       if (sekunde==60)                    %Nece biti tolika greska kao sto bi bila za sate
          sekunde=0;                   %Ako je 60 sekundi,napisi da je to nula
       end
end

sati=fix(UgaoSata/30);   %12 mogucnosti za sate,znaci delimo ugao sa 30
if (sati==0)             %Ako je nula sati,napisi da je 12 sati
    sati=12;
end
minuti=round(UgaoMin/6); %60 mogucnosti za minute,znaci delimo sa 6
if (minuti==60)          %Ako je 60,napisi da je nula
    minuti=0;
end

%Razlike u kodu nadogradjene funkcije
%Odavde samo za plotovanje

if (length(linije)==1)
    xy1=[linije(1).point1;linije(1).point2];
    figure(); imshow(I); hold on;
    plot([xy1(1,1) xy1(2,1)],[xy1(1,2) xy1(2,2)],'Color','y','LineWidth', 5); hold off;
elseif (length(linije)==2)
    xy1=[linije(1).point1;linije(1).point2];
    xy2=[linije(2).point1;linije(2).point2];
    figure(); imshow(I); hold on;
    plot([xy1(1,1) xy1(2,1)],[xy1(1,2) xy1(2,2)],'Color','y','LineWidth', 5); hold on;
    plot([xy2(1,1) xy2(2,1)],[xy2(1,2) xy2(2,2)],'Color','y','LineWidth', 5); hold off;
end
if (length(linije)==3)
    xy1=[linije(1).point1;linije(1).point2];
    xy2=[linije(2).point1;linije(2).point2];
    xy3=[linije(3).point1;linije(3).point2];
    figure(); imshow(I); hold on;
    plot([xy1(1,1) xy1(2,1)],[xy1(1,2) xy1(2,2)],'Color','y','LineWidth', 5); hold on;
    plot([xy2(1,1) xy2(2,1)],[xy2(1,2) xy2(2,2)],'Color','y','LineWidth', 5); hold on;
    plot([xy3(1,1) xy3(2,1)],[xy3(1,2) xy3(2,2)],'Color','y','LineWidth', 5); hold off;
end
    
end