
function [PlaveBonus,CrveneBonus]=extract_dice_score_bonus(I)

%
%Ime funkcije:extract_dice_score_bonus
%
%Funkcija prima kao ulazni argument sliku kockica.Na slici se moze nalaziti
%proizvoljan broj kockica na kojima se mogu nalaziti brojevi prikazani
%plavim i crvenim kruzicima.Ova funkcija vraca informaciju o tome koliko je
%broj bacen na kockicama sa plavim brojevima i koliki je broj bacen na
%kockicama sa crvenim brojevima,odnosno kruzicima,pri cemu vraca informaciju
%o tome koliko je zapravo kockica baceno i koliki je tacan broj bacen na tim
%kockicama. Ova funkcija koristi to sto je prosledjena slika u RGB formatu,
%pa na osnovu informacija o plavoj i crvenoj boji nalazi odgovarajuci prag
%i preko njega detektuje plave i crvene kruzice.Nakon detektovanja praga,
%nalazimo krugove i njihove centre.Prostim prebrojavanjem broja centara plavih i crvenih krugova
%dobijamo informaciju o bacenom broju na kockicama.
%
%Izgled funkcije:
%
%        [PlaveBonus, CrveneBonus]=extract_dice_score_bonus(I)
%
%Funkcija prima ulaznu sliku I u RGB formatu,a vraca nizove,odnosno informacije o broju
%plavih i crvenih kruzica bacenih na svakoj kockici,odnosno bacenih brojeva na kockicama.
%
%Funkcija nema podrazumevanih vrednosti,vec joj se eksplicitno prosledjuje
%slika u RGB formatu,a kao izlaz moze vratiti prazne nizove ako na
%trenutnoj slici nema uopste plavih ili crvenih brojeva,odnosno kruzica.
%
%Primer:
%
%          I=imread('dices1.jpg');
%          [PlaveBonus,CrveneBonus]=extract_dice_score_bonus(I);
%
%
%See also: extract_dice_score
%
%Funkcija extract_dice_score_bonus predstavlja nadogradnju vec postojece
%funkcije exctract_dice_score koja vraca nizove za plave i crvene kockice i
%vraca informaciju o tome koliko je plavih,a koliko crvenih kockica baceno 
%i koji se broj na njima nalazi
%
%
% Dan kreacije: 28.12.2019. (Petkovic Uros)
% Poslednje izmene: 28.12.2019. (Petkovic Uros)
%


%Prvo radimo obradu za plave kockice i delimo sliku u 3 plejna
R=I(:,:,1);
G=I(:,:,2);
B=I(:,:,3);   %Iz tehnickih razloga prolazimo kroz veci filtar ovog puta kako ne bismo uzeli u obzir
R=medfilt2(R,[9 9]); %neke segmente slike koji nam nisu od znacaja,npr. brojevi koji se naziru na 
%bocnim stranama kocke koji ne treba da se ocitaju
Plave=1.*(R<51 & G<75 & B>70);  %Postavljamo 3 praga eksperimentalno uvidom u originalnu sliku uz
%pomoc naredbe Data Cursor
[center1, radius] = imfindcircles(Plave, [2, 150]);  %Nalazimo centre kruzica
center11 = find(radius>3);  %Uzimamo one ciji je radijus veci od 3
figure(1); imshow(I);  %Prikazujemo sliku sa plavim kruzicima
viscircles(center1(center11,:), radius(center11), 'Color', 'b', 'LineWidth', 3);
hold on;

%Za crvene kockice radimo istu obradu kao i prethodni put,izdvajamo sva
%tri plejna
R=I(:,:,1);
G=I(:,:,2);
B=I(:,:,3);
R=medfilt2(R,[7 7]);  %Prolazimo kroz filtre radi dodatne obrade
G=medfilt2(G,[7 7]);
B=medfilt2(B,[7 7]);
Crvene=1.*(R>100 & G<90 & B<106);     %Postavljamo prag
[center2, radius] = imfindcircles(Crvene, [2, 150]);  %Nalazimo centre krugova
center22 = find(radius>3);  %Uzimamo one ciji je radijus veci od 3
figure(1);        %Prikazujemo ukupnu sliku sa plavim i crvenim kruzicima
viscircles(center2(center22,:), radius(center22), 'Color', 'r', 'LineWidth', 3);
%U ovom delu koda u stvari nalazimo koliko je bilo baceno plavih,a koliko
%crvenih kockica i koliki je broj na njima bio bacen.Najpre pravimo prazne
%nizove za Plave Bonus i za flag duzine jednake broju plavih kruzica.Ideja
%je prolaziti kroz dva fora i porediti kruzice,ako se oni nalaze unutar
%nekog praga,odnosno dovoljno su blizu,onda datu sekciju povecaj za
%1,odnosno prvi put povecaj za 2 kako bi uracunao i prvi kruzic pri cemu
%postavljamo dati flag na 1 ako smo taj kruzic vec uparili sa nekim od
%kruzica.Ako je taj flag postavljen na 1,onda u sledecoj iteraciji uopste
%necemo uzimati taj kruzic u razmatranje,vec cemo ga odmah preskociti.Ovde
%se javlja i ogranicen slucaj kada je na kockici samo 1 broj i nema oko
%sebe blizih brojeva,onda se pita da li je data sekcija prazna i da li je
%flag nula,odnosno nisi ga dodelio nikome,onda tu sekciju povecaj za 1.Isti
%kod ponavljan je i za crvene kockice.Kao ideja koriscen je klasican kod za
%sortiranje naucen na predmetu Algoritmi i strukture podataka.
PlaveBonus=zeros(1,numel(center11));   %Prazni nizovi
flag=zeros(1,numel(center11));
for i=1:numel(center11)-1   %Od prvo do pretposlednjeg idem
    if(flag(i)==0)          %Od drugog,odnosno trenutnog do poslednjeg idem
        for k=i+1:numel(center11)   %Idem kroz sve kruzice i proveravam im razdaljine
            if(sqrt(((center1(k,1)-center1(i,1))^2+(center1(k,2)-center1(i,2))^2))<51)
                if (PlaveBonus(i)==0)  %Ako su blizu i do tad nisam imao jos nikoga,povecaj za 2
                    PlaveBonus(i)=PlaveBonus(i)+2;
                    flag(k)=1;
                else                  %U suprotnom povecaj za 1
                    PlaveBonus(i)=PlaveBonus(i)+1; 
                    flag(k)=1;        %U oba slucaja postavi flag na 1 kako ne bi ponovo uzeo 
                end                   %dati kruzic u obzir
            end
        end
    end
end
for i=1:numel(center11)       %Istureni slucaj kada imamo samo jedan broj na kockici
    if (PlaveBonus(i)==0 && flag(i)==0)  %Prazna sekcija i flag na 0,znaci imamo jedan broj
    PlaveBonus(i)=PlaveBonus(i)+1; %Povecaj tu sekciju za 1
    end
end
PlaveBonus=PlaveBonus(PlaveBonus>0);  %Posto imamo niz koji je duzine onoliko koliko smo imali
%plavih kruzica,hocu da mi funkcija vrati samo one sekcije na kojima se
%nalaze kruzici,a ostale nule necu da mi posalje
CrveneBonus=zeros(1,numel(center22));  %Isti postupak ponavljam i za crvene
flag1=zeros(1,numel(center22));        %Pravim nizove
for i=1:numel(center22)-1              %Idem kroz dva fora
    if(flag1(i)==0)
        for k=i+1:numel(center22)      %Ako su u opsegu,odnosno blizu su,ulaze u if
            if(sqrt(((center2(k,1)-center2(i,1))^2+(center2(k,2)-center2(i,2))^2))<50)
                if (CrveneBonus(i)==0) %Ako do sad nisi nasao nijednog,povecaj za 2
                     CrveneBonus(i)=CrveneBonus(i)+2;
                    flag1(k)=1;    %Postavi flag
                else
                    CrveneBonus(i)=CrveneBonus(i)+1;  %Inace povecaj za 1
                    flag1(k)=1;             %Svakako postavi flag
                end
            end
        end
    end
end
for i=1:numel(center22)   %Istureni slucaj sa kockicom sa jednim brojem
    if (CrveneBonus(i)==0 && flag1(i)==0) 
    CrveneBonus(i)=CrveneBonus(i)+1;  %Povecaj za jedan
    end
end
CrveneBonus=CrveneBonus(CrveneBonus>0); %Hocu da mi vratis samo one nenulte sekcije
%Ovo su izlazi funkcije,odnosno nizovi koji nam govore o tome koliko je
%bilo bacenih kockica i koliki je broj na njima detektovan
end
