
function J=non_local_means(I,K,kij,var,h)

%Ime funkcije: non_local_means
%
%Funkcija se koristi za nelokalno usrednjavanje slike. Funkcija ima sledeci
%izgled:
%               J=non_local_means(I,K,S,var,h)
%
% gde je prvi parametar ulazna slika I u double formatu, drugi parametar K
% dimenzija susedstva piksela,treci parametar S dimenzija oblasti
% pretrazivanja piksela,cetvrti argument procenjena varijansa suma u
% zadatoj ulaznoj slici i poslednji,peti argument,h faktor odsumljivanja.
%Slika u svakom trenutku mora imati prosledjen svaki od argumenata,inace
%je poziv funkcije nepotpun i ona nije u stanju da se izvrsi i vraca
%odgovarajucu gresku,sto je ekvivalentno tome da funkcija nema
%podrazumevane vrednosti.
%
%Opis funkcije: Ova funkcija namenjena je nelokalnom usrednjavanju slike,
%tacnije,za odsumljavanje zasumljene slike poznate varijanse.
%Naime,funkcija dobija odredjene ulazne parametre,dimenzije prozora
%susedstva i prozora oblasti pretrazivanja,kao i varijansu suma i stepen
%odsumljavanja. Od stepena odsumljavanja zavisi koliko ce u stvari biti
%dobro otklonjen sum,pored ostalih parametara.U svakoj oblasti
%pretrazivanja SxS,koja je veca od susedstva KxK,vrsi se obrada slike tako
%sto se ide po celoj oblasti pretrazivanja tako da se obuhvati svaki piksel
%u datoj oblasti prozorom susedstva nalazeci se u sredini njega.Za svaki
%piksel racunamo njegovo euklidsko rastojanje u odnosu na ostale piksele
%obuhvacene oblasti pretrazivanja i taj medjurezultat koristimo u daljem
%toku programa tako sto on figurise u promenljivoj "wij" koju mnozimo sa
%trenutnim pikselom u oblasti pretrazivanja i smestamo u promenljivu NL,na
%kraju prolaska kroz datu oblast pretrazivanja vrsimo normalizaciju
%promenljive NL deleci je sumom svih elemenata sij,odnosno sumom S,cime smo
%dobili novu vrednost jednog od piksela.Ova pravilnost se ponavlja kroz
%celu sliku dok ne prodjemo kroz sve piksele i ne dobijemo novu,izlaznu
%sliku.
%
% Primer:
%
% I=im2double(imread('hugo.png'));
% figure; imshow(I);
% set(gcf, 'Name', 'Ulazna slika');
% J=non_local_means(I,3,33,0.0075,0.05);
% figure; imshow(J);
% set(gcf, 'Name', 'Izlazna slika');
%
%Sama obrada je zametna,zahteva dosta vremena,pogotovu kada su u pitanju
%prozori vecih dimenzija,stoga treba imati strpljenja.Ovo,takodje,zavisi i
%od brzine uredjaja na kojem se vrsi data obrada.
%
%
% See also: imnmfilt
%
% Dan kreacije: 29.11.2019. (Petkovic Uros)
% Poslednje izmene: 29.11.2019. (Petkovic Uros)


if (nargin ~=5)                             %Ako nismo prosledili tacan broj argumenata
   disp('Pogresan broj parametara');        %prijavljujemo odgovarajucu gresku
   J=0;
   return
end

if (iseven(K)) || (iseven(kij))                 %Ako unesemo pogresnu vrednost za K i S
   disp('K i S nisu neparni brojevi');        %prijavljujemo odgovarajucu gresku
   J=0;
   return
end

if (isa(I,'double')==0)                     %Ako nije prosledjena ulazna slika u double formatu
    disp('Slika nije u double formi');      %funkcija vraca odgovarajucu gresku
    J=0;
    return
end

[M,N] = size(I);        %Dimenzije ulazne slike
f=(K-1)/2;              %Prozor susedstva (2f+1)*(2f+1)=KxK
t=(kij-1)/2;              %Prozor oblasti pretrazivanja  (2t+1)*(2t+1)=SxS
 

J = zeros(M,N);         %Pravimo praznu matricu za smestanje izlazne slike
I = padarray(I,[f,f],'symmetric');  %Prosirujemo ulaznu sliku u ogledalu za vrednosti f zato sto
                                    %je to ceo deo polovine dimenzije prozora K,tako da bi se prvi
                                    %piksel u okviru prve oblasti
                                    %pretrazivanja mogao naci u sredini
                                    %datog prozora susedstva
%Ovo smo mogli uraditi i na drugi nacin prosirivanjem ulazne slike za vrednosti
%t,odnosno ceo deo polovine dimenzije S,pri cemu prosirujemo sliku za mnogo
%vise piksela,subjektivno je izabran prvi slucaj jer smatram da su
%rezultati dovoljno dobri i za prvi izbor,a zbog smanjenih dimenzija slike
%izvrsavanje koda koji je vec sam po sebi dovoljno dug znatno smanjeno

K2=K*K;    %Broj elemenata u prozoru susedstva koji koristimo u formuli za
%promenljivu sij kao koeficijent ispred eksponencijalnog clana kao u
%formuli sa predavanja

for i=1:M           %Idemo po svim poljima matrice ulazne slike
    for j=1:N       %Imamo pristup pikselu i,j
        im = i+f;   %Zbog produzenih dimenzija ulazne slike i centriranja zadatog piksela
        jn= j+f;    %u centar kvadratnog prozora pomeramo trenutnu vrednost piksela za f
        
                    %B1 prozor susedstva prvog piksela u odnosu na koji
                    %vrsimo usrednjavanje
        B1 = I(im-f:im+f , jn-f:jn+f);  %Kvadratno susedstvo KxK
        
                    %Odredjujemo pocetne i krajnje granice oblasti
                    %pretrazivanja u odnosu na posmatrani piksel
        poc = max(im-t, f+1);
        kraj = min(im+t, M+f);
        poc1 = max(jn-t, f+1);
        kraj1 = min(jn+t, N+f);
                    %Inicijalizujemo vrednosti 
        Ikapa=0;                        %Procena datog piksela pre normalizacije
        kij =0;                         % Suma svih wij pojacanja
                                        % For koji prolazi kroz sve piksele u oblasti pretrazivanja
        for r=poc:kraj                  %i racuna euklidsko rastojanje u odnosu na dati pocetni piksel
            for s=poc1:kraj1            %a potom racuna pojacanje wij kojim se mnozi trenutni piksel
                                        %koji se kasnije normalizuje
                                        %deljenjem sa kij
                                        %B2 kvadratno susedstvo trenutnog
                                        %piksela sa kojim uporedjujemo
                                        %glavni piksel oblasti
                                        %pretrazivanja
                                        
                B2 = I(r-f:r+f, s-f:s+f);
                                        %Suma pojacane euklidske distance
                                        %posmatranog i trenutnog piksela u
                                        %oblasti pretrazivanja pomnozene
                                        %reciprocnom vrednoscu broja
                                        %elemenata u kvadratnom susedstvu
                d2 =(1/K2)*sum(sum((B1-B2).*(B1-B2)));
                                        %Nakon toga racunamo wij po datoj
                                        %formuli
                wij = exp(-max(d2-2*var,0)./(h*h));
                                        %Apdejtujemo kij i Ikapa
                kij = kij + wij;
                Ikapa = Ikapa + wij*I(r,s);
            end
        end
                                        %Normalizujemo Ikapa i smestamo
                                        %dobijenu novu vrednost piksela u
                                        %izlaznu matricu i ovaj postupak
                                        %ponavljamo za svaki piksel u
                                        %ulaznoj slici
        J(i,j) = Ikapa/kij;             %Krajnji ishod je odsumljena slika
    end
end
end