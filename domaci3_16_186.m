
%% Prvi zadatak

clc; clear all; close all;

%Za svaku datu sliku pokrecemo funkcije extract_time i extract_time_bonus
%Funkcija extract_time vraca dve vrednosti,sate i minute,a onda mi te
%rezultate prikazujemo na komandnoj liniji u zavisnosti od sata koji se
%nalazi na datoj slici.Funkcija extract_time_bonus nam vraca tri
%vrednosti,sate,minute i sekunde,s tim sto ce promenljiva sekunde biti
%jednaka -1 ako smo prosledili sliku na kojoj ne postoji kazaljka za
%sekunde,a ako ona postoji,onda vraca odgovarajuci broj sekundi.Siri opis
%datih funkcija bice opisan u okviru samih zaglavlja funkcija. U ovom kodu
%pokretali smo funkciju za svih 12 slika i prikazivali rezultate.Potrebno
%je napomenuti da dati parametri u ovim funkcijama odgovaraju mom
%uredjaju,odnosno kompjuteru,i da se(znajuci iz iskustva jer sam domaci
%radio i na fakultetu na kompjuteru i kod kuce) moze desiti da ne radi za
%sve slike na nekom drugom uredjaju,a kao uverenje u to da su funkcije
%radile za sve primere bice date slike u izvestaju uz detaljan opis
%postupka.

%Prvi sat
I1=im2double(rgb2gray(imread('clock1.png')));
[sati1,minuti1]=extract_time(I1);
disp(['Prvi sat: ', num2str(sati1), ' :  ' ,num2str(minuti1)]);
[sati_1,minuti_1,sekunde_1]=extract_time_bonus(I1);
disp(['Prvi sat bonus: ', num2str(sati_1), ' :  ' ,num2str(minuti_1),' : ',num2str(sekunde_1)]);
%Drugi sat
I2=im2double(rgb2gray(imread('clock2.png')));
[sati2,minuti2]=extract_time(I2);
disp(['Drugi sat: ', num2str(sati2), ' :  ' ,num2str(minuti2)]);
[sati_2,minuti_2,sekunde_2]=extract_time_bonus(I2);
disp(['Drugi sat bonus: ', num2str(sati_2), ' :  ' ,num2str(minuti_2),' : ',num2str(sekunde_2)]);
%Treci sat
I3=im2double(rgb2gray(imread('clock3.png')));
[sati3,minuti3]=extract_time(I3);
disp(['Treci sat: ', num2str(sati3), ' :  ' ,num2str(minuti3)]);
[sati_3,minuti_3,sekunde_3]=extract_time_bonus(I3);
disp(['Treci sat bonus: ', num2str(sati_3), ' :  ' ,num2str(minuti_3),' : ',num2str(sekunde_3)]);
%Cetvrti sat
I4=im2double(rgb2gray(imread('clock4.png')));
[sati4,minuti4]=extract_time(I4);
disp(['Cetvrti sat: ', num2str(sati4), ' :  ' ,num2str(minuti4)]);
[sati_4,minuti_4,sekunde_4]=extract_time_bonus(I4);
disp(['Cetvrti sat bonus: ', num2str(sati_4), ' :  ' ,num2str(minuti_4),' : ',num2str(sekunde_4)]);
%Peti sat
I5=im2double(rgb2gray(imread('clock5.png')));
[sati5,minuti5]=extract_time(I5);
disp(['Peti sat: ', num2str(sati5), ' :  ' ,num2str(minuti5)]);
[sati_5,minuti_5,sekunde_5]=extract_time_bonus(I5);
disp(['Peti sat bonus: ', num2str(sati_5), ' :  ' ,num2str(minuti_5),' : ',num2str(sekunde_5)]);
%Sesti sat
I6=im2double(rgb2gray(imread('clock6.png')));
[sati6,minuti6]=extract_time(I6);
disp(['Sesti sat: ', num2str(sati6), ' :  ' ,num2str(minuti6)]);
[sati_6,minuti_6,sekunde_6]=extract_time_bonus(I6);
disp(['Sesti sat bonus: ', num2str(sati_6), ' :  ' ,num2str(minuti_6),' : ',num2str(sekunde_6)]);
%Sedmi sat
I7=im2double(rgb2gray(imread('clock7.jpg')));
[sati7,minuti7]=extract_time(I7);
disp(['Sedmi sat: ', num2str(sati7), ' :  ' ,num2str(minuti7)]);
[sati_7,minuti_7,sekunde_7]=extract_time_bonus(I7);
disp(['Sedmi sat bonus: ', num2str(sati_7), ' :  ' ,num2str(minuti_7),' : ',num2str(sekunde_7)]);
%Osmi sat
I8=im2double(rgb2gray(imread('clock8.jpg')));
[sati8,minuti8]=extract_time(I8);
disp(['Osmi sat: ', num2str(sati8), ' :  ' ,num2str(minuti8)]);
[sati_8,minuti_8,sekunde_8]=extract_time_bonus(I8);
disp(['Osmi sat bonus: ', num2str(sati_8), ' :  ' ,num2str(minuti_8),' : ',num2str(sekunde_8)]);
%Deveti sat
I9=im2double(rgb2gray(imread('clock9.jpg')));
[sati9,minuti9]=extract_time(I9);
disp(['Deveti sat: ', num2str(sati9), ' :  ' ,num2str(minuti9)]);
[sati_9,minuti_9,sekunde_9]=extract_time_bonus(I9);
disp(['Deveti sat bonus: ', num2str(sati_9), ' :  ' ,num2str(minuti_9),' : ',num2str(sekunde_9)]);
%Deseti sat
I10=im2double(rgb2gray(imread('clock10.png')));
[sati10,minuti10]=extract_time(I10);
disp(['Deseti sat: ', num2str(sati10), ' :  ' ,num2str(minuti10)]);
[sati_10,minuti_10,sekunde_10]=extract_time_bonus(I10);
disp(['Deseti sat bonus: ', num2str(sati_10), ' :  ' ,num2str(minuti_10),' : ',num2str(sekunde_10)]);
%Jedanaesti sat
I11=im2double(rgb2gray(imread('clock11.png')));
I11=I11(23:282,25:280);
[sati11,minuti11]=extract_time(I11);
disp(['Jedanaesti sat: ', num2str(sati11), ' :  ' ,num2str(minuti11)]);
[sati_11,minuti_11,sekunde_11]=extract_time_bonus(I11);
disp(['Jedanaesti sat bonus: ', num2str(sati_11), ' :  ' ,num2str(minuti_11),' : ',num2str(sekunde_11)]);
%Dvanaesti sat
I12=im2double(rgb2gray(imread('clock12.jpg')));
[sati12,minuti12]=extract_time(I12);
disp(['Dvanaesti sat: ', num2str(sati12), ' :  ' ,num2str(minuti12)]);
[sati_12,minuti_12,sekunde_12]=extract_time_bonus(I12);
disp(['Dvanaesti sat bonus: ', num2str(sati_12), ' :  ' ,num2str(minuti_12),' : ',num2str(sekunde_12)]);

%% Drugi zadatak

clc; clear all; close all;

%U ovom kodu pokrecemo dve funkcije extract_dice_score i
%extract_dice_score_bonus za 12 datih primera kockica.Prva funkcija treba
%da vrati broj plavih i crvenih kruzica na kockicama i da izvrsi njihov
%zbir na kraju kao uvid(ovo je subjektivno resenje).Druga funkcija treba da
%vrati niz brojeva za plave i za crvene kockice takav da pokaze koliko je
%bilo tacno plavih i crvenih kruzica na kojoj kockici.Ovaj program radi za
%sve vrednosti i za bilo koji broj kockica.Ako imamo npr. 15 zbir,teorijski
%je moguce da imamo 15 kockica,sto je ovim kodom ispostovano.Sve slike sa
%tacno oznacenim krugovima plave i crvene boje,kao i rezultati funkcija
%bice prikazani u izvestaju,a detaljan opis funkcija u samim zaglavljima
%funkcija i,takodje,u izvestaju.Jedna od napomena je kao sto je to bilo i
%za prvi zadatak da ovi rezultati rade za sve primere na mom uredjaju i da
%se mozda moze desiti da na nekom drugom ne bude isto posto sam zadatak resavao
%preko RGB komponenti,a ovaj kolor sistem je zavisan od uredjaja,ali nadam
%se da to ne bi trebalo da bude slucaj.
%Prva slika kockica
I1=imread('dices1.jpg');
[Plave1,Crvene1]=extract_dice_score(I1);
Zbir1=Plave1+Crvene1;
[PlaveBonus1,CrveneBonus1]=extract_dice_score_bonus(I1);
%Druga slika kockica
I2=imread('dices2.jpg');
[Plave2,Crvene2]=extract_dice_score(I2);
Zbir2=Plave2+Crvene2;
[PlaveBonus2,CrveneBonus2]=extract_dice_score_bonus(I2);
%Treca slika kockica
I3=imread('dices3.jpg');
[Plave3,Crvene3]=extract_dice_score(I3);
Zbir3=Plave3+Crvene3;
[PlaveBonus3,CrveneBonus3]=extract_dice_score_bonus(I3);
%Cetvrta slika kockica
I4=imread('dices4.jpg');
[Plave4,Crvene4]=extract_dice_score(I4);
Zbir4=Plave4+Crvene4;
[PlaveBonus4,CrveneBonus4]=extract_dice_score_bonus(I4);
%Peta slika kockica
I5=imread('dices5.jpg');
[Plave5,Crvene5]=extract_dice_score(I5);
Zbir5=Plave5+Crvene5;
[PlaveBonus5,CrveneBonus5]=extract_dice_score_bonus(I5);
%Sesta slika kockica
I6=imread('dices6.jpg');
[Plave6,Crvene6]=extract_dice_score(I6);
Zbir6=Plave6+Crvene6;
[PlaveBonus6,CrveneBonus6]=extract_dice_score_bonus(I6);
%Sedma slika kockica
I7=imread('dices7.jpg');
[Plave7,Crvene7]=extract_dice_score(I7);
Zbir7=Plave7+Crvene7;
[PlaveBonus7,CrveneBonus7]=extract_dice_score_bonus(I7);
%Osma slika kockica
I8=imread('dices8.jpg');
[Plave8,Crvene8]=extract_dice_score(I8);
Zbir8=Plave8+Crvene8;
[PlaveBonus8,CrveneBonus8]=extract_dice_score_bonus(I8);
%Deveta slika kockica
I9=imread('dices9.jpg');
[Plave9,Crvene9]=extract_dice_score(I9);
Zbir9=Plave9+Crvene9;
[PlaveBonus9,CrveneBonus9]=extract_dice_score_bonus(I9);
%Deseta slika kockica
I10=imread('dices10.jpg');
[Plave10,Crvene10]=extract_dice_score(I10);
Zbir10=Plave10+Crvene10;
[PlaveBonus10,CrveneBonus10]=extract_dice_score_bonus(I10);
%Jedanaesta slika kockica
I11=imread('dices11.jpg');
[Plave11,Crvene11]=extract_dice_score(I11);
Zbir11=Plave11+Crvene11;
[PlaveBonus11,CrveneBonus11]=extract_dice_score_bonus(I11);
%Dvanaesta slika kockica
I12=imread('dices12.jpg');
[Plave12,Crvene12]=extract_dice_score(I12);
Zbir12=Plave12+Crvene12;
[PlaveBonus12,CrveneBonus12]=extract_dice_score_bonus(I12);


%% Treci zadatak

clc; clear all;close all;
%U ovom kodu pokrecemo funkciju canny_edge_detection za slike
%lena.tif,van.tif,camerman.tif i house.tif.Ideja ove funkcije je izdvajanje
%ivica Kanijevim algoritmom.Najpre smo filtrirali sliku Gausovim filtrom,
%odredili magnitudu i ugao gradijenta,kvantizovali ugao na 4
%vrednosti,potisnuli lokalne nemaksimume,odredili jake i slabe ivice i
%nakon toga u jake ivice ukljucili i one slabe ivice koje su povezane sa
%nekim jakim ivicama.Rezultate funkcije za date parametre bice prikazani u
%izvestaju.
I1 = im2double(imread('lena.tif'));
i=1;
for sigma1=[0.834 2 2.8]
  for low1=[0.02 0.05 0.1]
    for high1=[0.05 0.1 0.23]
       J1=canny_edge_detection(I1,sigma1,low1,high1);
       figure(i); imshow(J1);
       title(['sigma1=',num2str(sigma1),'low1',num2str(low1),'high1',num2str(high1)]);
       i=i+1;
    end
  end
end

I2 = im2double(imread('camerman.tif'));
sigma2=1.5;
low2=0.06;
high2=0.12;
J2=canny_edge_detection(I2,sigma2,low2,high2);
figure(); imshow(J2);

I3 = im2double(imread('van.tif'));
sigma3=1.4;
low3=0.046;
high3=0.058;
J3=canny_edge_detection(I3,sigma3,low3,high3);
figure(); imshow(J3);

I4 = im2double(imread('house.tif'));
sigma4=1.2;
low4=0.017;
high4=0.065;
J4=canny_edge_detection(I4,sigma4,low4,high4);
figure(); imshow(J4);