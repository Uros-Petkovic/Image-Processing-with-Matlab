function  J = canny_edge_detection(I,sigma,low,high)
%   
%Ime funkcije
%
%Funkcija kao argumente prima ulaznu sliku,standardnu devijaciju koju treba upotrebiti za
%filtriranje slike Gausovim filtrom,a,takodje,prima i apsolutne vrednosti nizeg
%i viseg praga koji nisu vezani za maksimalni gradijent u slici.Ova
%funkcija ima za ulogu da nadje ivice na slici koriscenjem Kanijevog
%algoritma.Najpre pocetnu sliku isfiltrira Gausovim filtrom.U okviru ove
%funkcije koriscena je funkcija gradient kojoj se prosledjuju informacije o
%Gausovom filtru,tako da samim tim nalazi odgovarajuce Hx i Hy parametre
%takve da u nastavku moze naci odgovarajucu magnitudu i ugao gradijenta
%filtrirane slike Gausovim filtrom.Nakon toga vrsi kvantizaciju ugla
%gradijenta na cetiri vrednosti,0,45,-45 i 90 stepeni.Potisuje lokalne
%nemaksimume,nalazi matrice slabih i jakih ivica.Kao konacno resenje uzima
%se u obzir ono resenje u koje su ukljucene i slabe ivice koje su povezane
%sa nekom od jakih ivica.Tako dobijena slika salje se na izlaz funkcije.Ona
%predstavlja mapu jakih ivica slike.
%
%Izgled funkcije:
%
%       J = canny_edge_detection(I,sigma,low,high);
%
% Gde je J izlazna slika jakih ivica,I ulazna slika,sigma standardna
% devijacija,low i high donji i gornji prag za gradijent.
%
%Funkcija nema podrazumevane vrednosti vec joj se eksplicitno moraju
%proslediti sva tri argumenta.
%
%Primer:
%
%           I2 = im2double(imread('camerman.tif'));
%           sigma2=1.5;
%           low2=0.06;
%           high2=0.12;
%           J2=canny_edge_detection(I2,sigma2,low2,high2);
%           figure(); imshow(J2);
%
%
%See also: canny,edge
%
%
% Dan kreacije: 28.12.2019. (Petkovic Uros)
% Poslednje izmene: 28.12.2019. (Petkovic Uros)
%
%


    %K parametar nalazim po uputstvu zadatka,odnosno informaciji da k mora
    %biti barem 6 puta veci od prosledjene standardne devijacije zaokruzen
    %na prvi sledeci ceo broj
    k=ceil(6*sigma);        %Zaokruzuje na veci ceo broj

    if (2*floor(k/2)==k) %Ako je paran broj povecaj ga za jedan
        k=k+1;           %Ako je broj paran,mora se povecati za jedan jer fspecial 
    end                  %prima samo neparne argumente
    if (k<3)             %Ako je k<3 uzimamo prozor 3x3 jer necemo 1x1 prozor,3x3 je najmanji
        k=3;
    end
    h=fspecial('gaussian',[k k], sigma);  %Gausov filtar pravimo
    [Hx,Hy]=gradient(h);       %Koriscenje Gausovog filtra za nalazenje Hx i Hy
    Hx = Hx/sum(abs(Hx(:)));   %Za vertikalni gradijent
    Hy = Hy/sum(abs(Hy(:)));   %Za horizontalni gradijent
    Gx = imfilter(I, Hx, 'replicate', 'same'); %Nalazimo Gx komponentu
    Gy = imfilter(I, Hy, 'replicate', 'same'); %Nalazimo Gy komponentu
    figure('Name', 'Vertikalni gradijent'); imshow(Gx, []);  %Prikazujemo ih
    figure('Name', 'Horizontalni gradijent'); imshow(Gy, []);
    Mag = (Gx.^2 + Gy.^2).^0.5;    %Magnituda gradijenta
    Mag (Mag == 0) = 0;            
    Ph = atan(Gy./Gx);             %Ugao gradijenta
    figure('Name', 'Magnituda gradijenta'); imshow(Mag, []); %Prikaz magnitude
    figure('Name', 'Ugao gradijenta'); imagesc(Ph);     %Prikaz ugla gradijenta
    [M,N]=size(Ph);
    
    %Kvantizacija gradijenta na jedan od 4 pravca,odnosno na nivoe
    %0,45,-45 i 90 stepeni,kako ne bismo dodatno komplikovali kod i povecavali
    %slozenost algoritma,iskoristicemo mogucnost Matlaba da generise izlaze
    %uz pomoc logickih operatora i tako nesto sto bi bilo u 4 if-a
    %spakujemo u jednu naredbu u zavisnosti od odredjenih pragova.Vrednost
    %ugla ce biti kvantizovana na 0 ako se nadje u opsegu od -22.5 do 22.5
    %stepeni,u 90 ako se nadje od -90 do -67.5 i 67.5 i 90 stepeni,u -45
    %ako je u opsegu od -67.5 do -22.5 i konacno u 45 ako je u opsegu od
    %22.5 do 67.5
    
    Ph=-pi/4.*(Ph>=-3*pi/8 & Ph<-pi/8)+0.*(Ph>=-pi/8 & Ph<pi/8)+...
    pi/4.*(Ph>=pi/8 & Ph<3*pi/8)+pi/2.*((Ph>=-pi/2 & Ph<-3*pi/8) | (Ph>=3*pi/8 & Ph<=pi/2));
    
    Ph=Ph*180/pi;        %Prebacivanje u stepene
    Gs = padarray(Mag, [1, 1], 'replicate'); %Prosirujemo matricu za ram debljine 1 piksel
    %zato sto idemo prozorom 3x3 i hocu da mi u centru prozora bude moj
    %piksel,kako ne bih izasao iz opsega slike,celu sliku cu prosiriti za
    %ram velicine 1 piksel
    %Matrica koju koristimo za potiskivanje lokalnih nemaksimuma
    for i = 2:M+1           %U zavisnosti od toga u kom se kvantizovanom nivou ugla nalazimo
        for j = 2:N+1       %datu lokalnu 3x3 matricu mnozimo odgovarajucom 3x3 matricom suprotnog ugla
                if (Ph(i-1,j-1)==0)
                    k_lok = Gs(i-1:i+1,j-1:j+1).*[0 0 0;1 1 1;0 0 0];
                end
                if (Ph(i-1,j-1)==45)
                    k_lok = Gs(i-1:i+1,j-1:j+1).*[1 0 0;0 1 0;0 0 1];
                end
                if (Ph(i-1,j-1)==90)
                    k_lok = Gs(i-1:i+1,j-1:j+1).*[0 1 0;0 1 0;0 1 0];
                end
                if (Ph(i-1,j-1)==-45)
                    k_lok = Gs(i-1:i+1,j-1:j+1).*[0 0 1;0 1 0;1 0 0];
                end
                
                if(k_lok(2,2)~=max(k_lok(:)))  %Ako sredisnji element lokalne matrice 3x3 nije maksimum
                   Mag(i-1,j-1) = 0;  %onda tu vrednost odbacujemo,odnosno postavljamo vrednost magnitude
                end                   %za dati piksel na 0 jer je to lokalni nemaksimum
        end
    end
    %Nakon toga mozemo naci matrice jakih i slabih ivica,za to mozemo
    %iskoristiti informacije o visokom i niskom pragu,pri cemu cemo reci da
    %u jake ivice spadaju sve one kojima je magnituda veca od visokog
    %praga,a slabe one koje se nalaze izmedju niskog i visokog praga,pri
    %cemu cemo sve ostale odbaciti
    figure('Name', 'Potisnuti lokalni nemaksimumi'); imshow(Gs, []);
    J_jake=zeros(M,N);
    J_jake(Mag >= high)=1;       %Vece od praga,onda su jake
    J_slabe=zeros(M,N);
    J_slabe(Mag >= low & Mag < high)=1;  %Izmedju low i high slabe su
    figure('Name', 'Mapa jakih ivica'); imshow(J_jake, []); %Prikaz
    figure('Name', 'Mapa slabih ivica'); imshow(J_slabe, []);
        
    J_jake_2=zeros(M+2,N+2);      %Zbog prozora prosirujemo za 2 mesta
    J_jake_2(2:M+1,2:N+1)=J_jake; %Da bismo ukljucili i one slabe ivice oko jakih
    J_slabe_2=zeros(M+2,N+2);
    J_slabe_2(2:M+1,2:N+1)=J_slabe;
    flag=1;
    while flag                     
        flag=0;            %Ulaze u obzir i slabe ivice koje oko sebe imaju neku jaku ivicu,
        for i=2:M+1        %odnosno povezane su sa njom tako da se ta slaba ivica moze smatrati kao jaka
            for j=2:N+1
                if (J_slabe_2(i,j)==1)   %Ako je data ivica slaba
                    K_strong=J_jake_2(i-1:i+1,j-1:j+1); %3x3 prozor jaih ivica uzimam
                    if (sum(K_strong(:))~=0)   %Ako suma jakih nije jednaka nuli,onda se slaba ivica
                        J_jake_2(i,j)=1;  %moze proglasiti za jaku ivicu,postavljamo vrednost jake na 1
                        J_slabe_2(i,j)=0; %a vrednost koja je bila upisana u slabu ivicu postavljamo na 0
                        flag=1;
                    end
                end
                              
            end
        end
    end
    %Na kraju dobijamo konacan izgled slike jakih ivica istih dimenzija
    J=J_jake_2(2:M+1,2:N+1); %Krajnji rezultat jake ivice
    figure('Name', 'Slika izlaznih ivica'); imshow(J, []); %Prikazujemo datu sliku
      
    
end


