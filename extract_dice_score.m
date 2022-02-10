

function [PlaveUk, CrveneUk]=extract_dice_score(I)
%
%Ime funkcije:extract_dice_score
%
%Funkcija prima kao ulazni argument sliku kockica.Na slici se moze nalaziti
%proizvoljan broj kockica na kojima se mogu nalaziti brojevi prikazani
%plavim i crvenim kruzicima.Ova funkcija vraca informaciju o tome koliko je
%broj bacen na kockici sa plavim brojevima i koliki je broj bacen na
%kockicama sa crvenim brojevima,odnosno kruzicima,ne vracajuci informaciju
%o tome koliko je zapravo kockica baceno. Ova funkcija koristi to sto je
%prosledjena slika u RGB formatu,pa na osnovu informacija o plavoj i
%crvenoj boji nalazi odgovarajuci prag i preko njega detektuje plave i
%crvene kruzice.Nakon detektovanja praga,nalazimo krugove i njihove
%centre.Prostim prebrojavanjem broja centara plavih i crvenih krugova
%dobijamo informaciju o bacenom broju na kockicama.
%
%Izgled funkcije:
%
%        [PlaveUk, CrveneUk]=extract_dice_score(I)
%
%Funkcija prima ulaznu sliku I u RGB formatu,a vraca informaciju o broju
%plavih i crvenih kruzica,odnosno bacenih brojeva na kockicama.
%
%Funkcija nema podrazumevanih vrednosti,vec joj se eksplicitno prosledjuje
%slika u RGB formatu,a kao izlaz moze vratiti praznu promenljivu ako na
%trenutnoj slici nema uopste plavih ili crvenih brojeva,odnosno kruzica.
%
%Primer:
%
%          I1=imread('dices1.jpg');
%          [Plave1,Crvene1]=extract_dice_score(I1);
%          Zbir1=Plave1+Crvene1;
%
% Gde zbir predstavlja ukupan broj koji je bacen na kockicama.
%
%See also: extract_dice_score_bonus
%
%Funkcija extract_dice_score_bonus predstavlja nadogradnju vec postojece
%funkcije koja vraca nizove za plave i crvene kockice i vraca informaciju o
%tome koliko je plavih,a koliko crvenih kockica baceno i koji se broj na
%njima nalazi
%
%
% Dan kreacije: 28.12.2019. (Petkovic Uros)
% Poslednje izmene: 28.12.2019. (Petkovic Uros)
%

R=I(:,:,1);   %Prvo vrsimo obradu za plave kockice,delimo na tri plejna datu sliku
G=I(:,:,2);
B=I(:,:,3);
R=medfilt2(R,[7 7]);   %Prolazimo medijan filtrom kroz R
Plave= R>40;           %Odredjujemo prag kojim pronalazimo plave tacke na slici
Plave=1-Plave;         %Hocemo da nam bude bela slika na mestima kruzica,a crna u ostatku
[center1, radius] = imfindcircles(Plave, [2, 150]); %Nalazimo koordinate centara plavih krugova
center11 = find(radius>3);    %Uzimamo one krugove ciji je radijus veci od 3
figure(1); imshow(I);         %Prikazujemo sliku sa datim kruzicima
viscircles(center1(center11,:), radius(center11), 'Color', 'b', 'LineWidth', 3);
hold on;

%Za crvene kockice radimo isto,delimo sliku na 3 plejna,zbog kompleksnijeg
%nalazenja crvenih delova,u ovoj obradi koristimo sva tri plejna,pragovi su
%nadjeni eksperimentalno koristeci originalnu sliku i Data Cursor komandu
R=I(:,:,1);
G=I(:,:,2);
B=I(:,:,3);
R=medfilt2(R,[7 7]);
G=medfilt2(G,[7 7]);    %Prolazimo median filtrom
B=medfilt2(B,[7 7]);
Crvene=1.*(R>100 & G<90 & B<106);   %Postavljamo prag i nalazimo crvene delove
[center2, radius] = imfindcircles(Crvene, [2, 150]); %Nalazimo centre krugova
center22 = find(radius>3);   %Uzimamo one ciji je radijus veci od 3
figure(1);            %Prikazujemo sliku sa crvenim kruzicima i prethodno nadjenim plavim
viscircles(center2(center22,:), radius(center22), 'Color', 'r', 'LineWidth', 3);

PlaveUk=numel(center11); %Nalazimo broj plavih kruzica preko broja elemenata centara
CrveneUk=numel(center22);%Takodje,isto to radimo i za crvene kruzice
                         %Ove vrednosti su izlazi funkcije
end