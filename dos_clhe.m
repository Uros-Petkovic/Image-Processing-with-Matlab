
function J=dos_clhe(I,bit_depth,limit)
% Ime funkcije: dos_clhe
% Funkcija se koristi za poboljsanje kontrasta slike uz koriscenje
% ekvalizacije histograma slike sa mogucnoscu limitiranja histograma
% odredjenom vrednoscu koja se moze zadati.Limit se koristi kako bi se
% izbeglo pojacanje suma koji se moze pojaviti u slici.Funkcija maksimalno
% prima 3 argumenta,a moze se pozvati i sa manje argumenata.Tipican poziv
% funkcije izgleda ovako:
%
% J=dos_clhe(I,bit_depth,limit) gde je sa I oznacena ulazna slika koju
% treba obraditi,sa bit_depth oznacen je broj bita na kojim se slika moze
% predstaviti,kao i broj bita na kojima se pravi i sam histogram slike.On
% definise koliko slika moze imati vrednosti od 0 do 1.Ulazna slika se
% salje u opsegu od 0 do 1,a i izlazna slika je,takodje,u opsegu od 0 do 1.
% Limit limitira maksimalnu vrednost histograma skaliranog na opseg od 0 do
% 1,samim tim limitira pojacanje kontrasta.
%
% Podrazumevane vrednosti:
%
% Ako se funkciji prosledi samo jedan argument,slika funkcije,tada vrednost
% bit_depth uzima vrednost 8,a limit vrednost 0.01
%
% bit_depth=8; limit=0.01;
%
%Ako se funkciji proslede 2 argumenta,slika i bit_depth,onda limit uzima
%podrazumevanu vrednost 0.01
%
% limit=0.01;
% Ako se funkciji proslede neodgovarajuce vrednosti argumenata,kao sto je
% slike koja nije u double formatu,bit_depth koji je manji od 0 ili limit
% koji je van opega 0 do 1,funkcija vraca odgovarajucu gresku.
%
% Primer:
%
% I=imread('hugo.png');
% figure; imshow(I);
% set(gcf, 'Name', 'Ulazna slika');
% J=dos_clhe(I,8,0.01);
% figure; imshow(J);
% set(gcf, 'Name', 'Izlazna slika');
%
% See also: histeq, adapthisteq, hist, imhist
%
% Dan kreacije: 15.11.2019. (Petkovic Uros)
% Poslednje izmene: 15.11.2019. (Petkovic Uros)

if (nargin > 3) || (nargin < 1)             %Ako smo prosledili vise od 3 ili manje od 1
   disp('Pogresan broj parametara');        %argumenta,funkcija vraca odgovarajucu gresku
   J=0;
   return
else
    if nargin==2                            %Ako je broj prosledjenih argumenata 2,tada
         limit=0.01;                        %limit uzima podrazumevanu vrednost 0.01
    else 
        if nargin==1                        %Ako je broj prosledjenih argumenata 2,tada
            bit_depth=8;                    %limit uzima podrazumevanu vrednost 0.01,a
            limit=0.01;                     %bit_depth vrednost 8
        end
    end
end
if (isa(I,'double')==0)                     %Ako nije prosledjena ulazna slika u double formatu
    disp('Slika nije u double formi');      %funkcija vraca odgovarajucu gresku
    J=0;
    return
end

if (limit<0) || (limit>1)                   %Ako limit nije u opsegu od 0 do 1,funkcija vraca gresku
    disp('Limit ima nedozvoljenu vrednost');
    J=0;
    return
end

if (bit_depth < 0) || (bit_depth-round(bit_depth)~=0) %Ako bit_depth nema odgovarajucu vrednost,funkcija
    disp('Dubina bita ima nedozvoljenu vrednost');    %vraca odgovarajucu gresku
    J=0;
    return
end

if numel(size(I))==3                        %U zavisnosti od dimenzija slike,da li je 2D ili 3D,ulazimo u
   Ihsv=rgb2hsv(I);                         %jednu od if grana
   Iv=Ihsv(:,:,3);                          %U slucaju 3D slike,vrsicemo ekv. histograma u HSV kolor sistemu
   Ihist=zeros(1,2^bit_depth);              %Za pocetak pravimo prazan niz za histogram duzine 2^bit_depth
   N=numel(Iv);                             %Broj elemenata V matrice
   Is=round(Iv.*(2^bit_depth-1));
   for i=1:N
       Ihist(Is(i)+1)=Ihist(Is(i)+1)+1;     %Ako se komponenta slike nalazi na datom odbirku histograma,povecaj
   end                                      %njegovu vrednost za 1,na taj nacin dobijamo koliko se kojih komponenti
   Ihistn=Ihist./N;                         %slike nalazi na datom odbirku histograma,nakon cega njega delimo sa 
   summ=Ihistn(1,:)-limit;                  %ukupnim brojem elemenata slike,cime dobijamo opseg od 0 do 1
   summ(summ<0)=0;                          
   suma=sum(summ);                          %Nalazimo kolika je povrsina gornjeg dela koji treba odseci limitom
   suma=suma./2^(bit_depth);                %i nalazimo srednju vrednost koju cemo dodati svakom odbirku odsecenog
   Ihistn(Ihistn>limit)=limit;              %histograma,takodje,odsecamo histogram na max vrednost=limit
   Ihistn(1,:)=Ihistn(1,:)+suma;            %Dodajemo vrednost na svaki odbirak histograma
   figure; bar(Ihistn);                            
   title('Normalizovani histogram slike');  %Prikazujemo histogram radi uvida u trenutno stanje
   cdf = cumsum(Ihistn);                    %Nalazimo kumulativni histogram
   if (bit_depth>8) && (bit_depth<=16)
       lut=uint16(round(cdf*(2^bit_depth-1)));  %U zavisnosti od bit_depth pravimo odredjeni ekvivalent u integer domenu
       Iv=uint16(round(Iv*(2^bit_depth-1)));
   else
       if (bit_depth > 16)
           lut=uint32(round(cdf*(2^bit_depth-1)));
           Iv=uint32(round(Iv*(2^bit_depth-1)));
       else
           lut=uint8(round(cdf*(2^bit_depth-1)));
           Iv=uint8(round(Iv*(2^bit_depth-1)));
       end
   end
   Iv=intlut(Iv, lut);                           %Izvrsavamo preslikavanje ulaznih komponenti slike u izlazne uz pomoc
   Iv=double(Iv)./(2^bit_depth);                 %kumulativne sume
   Ihsv(:,:,3)=Iv;                               %Dodeljujemo promenjene komponente pocetnoj slici
   J=hsv2rgb(Ihsv);                              %Vracamo sliku u RGB format i saljemo je na izlaz
else
   Ihist=zeros(1,2^bit_depth);                   %Ako je slika u 2D formatu,sprovodimo isti postupak
   N=numel(I);                                   %samo sada obradjujemo crno-belu sliku
   Is=round(I.*(2^bit_depth));                   %Pravimo histogram,nalazimo koliku koji odbirak ima vrednost
   for i=1:N                               
     Ihist(Is(i)+1)=Ihist(Is(i)+1)+1;
   end
   Ihistn=Ihist./numel(I);                       %Delimo ga sa ukupnim brojem elemenata,nalazimo povrsinu odsecenog dela
   summ=Ihistn(1,:)-limit;                       %Odsecamo limit i dodajemo srednju vrednost svakom odbirku histograma
   summ(summ<0)=0;
   suma=sum(summ);
   suma=suma./2^(bit_depth);
   Ihistn(Ihistn>limit)=limit;
   Ihistn(1,:)=Ihistn(1,:)+suma;
   figure; bar(Ihistn);                            
   title('Normalizovani histogram slike');      %Prikazujemo histogram radi uvida u trenutno stanje
   cdf = cumsum(Ihistn);                         %Nalazimo kumulativni histogram
   if (bit_depth>8) && (bit_depth<=16)
     lut=uint16(round(cdf*(2^bit_depth-1)));     %Nalazimo odgovarajuci integer ekvivalent
     I=uint16(round(I*(2^bit_depth-1)));
   else
       if (bit_depth > 16)
        lut=uint32(round(cdf*(2^bit_depth-1)));
        I=uint32(round(I*(2^bit_depth-1)));
       else
        lut=uint8(round(cdf*(2^bit_depth-1)));
        I=uint8(round(I*(2^bit_depth-1)));
     end
   end
J1=intlut(I, lut);                               %Vrsimo preslikavanje ulaznih u izlazne komponente
J=im2double(J1);                                 %Saljemo izlaznu sliku u double formatu u opsegu od
end                                              %0 do 1 
  
end