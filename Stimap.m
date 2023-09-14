
function [p,c]=Stimap(xvect,nit)

%__________________________________________________________________________
%
% Function: Stimap
% ________________
%
%
%    Chiamata:  [p,c]=Stimap(xvect,nit);
%
%
%    Commento:
%
%    Stima ordine e costante asintotica di convergenza di un metodo 
%    iterativo per il calcolo degli zeri di una funzione utilizzando
%    le seguenti formule :
%           
%            | x_(k+1) - x_k |
%         ln ------------------
%            | x_k - x_(k-1) |                | x_(k+1) - x_k |
%    p = ---------------------------    c = -----------------------
%              | x_k - x_(k-1) |             | x_(k) - x_(k-1) |^p  
%        ln -----------------------
%            | x_(k-1) - x_(k-2) |
%
%
%
% Parametri in ingresso:
% ______________________
%
%     xvect      Vettore colonna di lunghezza nit, contenente 
%                le iterate ottenute da un metodo iterativo  
%                compresi gli eventuali punti di innesco.
%     nit        Lunghezza del vettore xvect.
%                
%
% Parametri in uscita:
% _____________________
%
%     p          Vettore colonna di lunghezza nit-1, contenente
%                le stime dell'ordine di convergenza calcolate 
%                ad ogni passo.
%     c          Vettore colonna di lunghezza nit-1, contenente 
%                le stime della costante asintotica di convergenza
%                calcolate ad ogni passo.
%__________________________________________________________________________



p([1,2],:)=zeros(2,1);
c([1,2],:)=zeros(2,1);

for i=3:nit-1
   diff1=abs(xvect(i+1)-xvect(i));	
   diff2=abs(xvect(i)-xvect(i-1));
   diff3=abs(xvect(i-1)-xvect(i-2));

   if abs(diff1)<=eps | abs(diff2)<=eps | abs(diff3)<=eps
      p(i)=p(i-1);
      c(i)=c(i-1);
   else
      num=log(abs(xvect(i+1)-xvect(i))/abs(xvect(i)-xvect(i-1)));
      den=log(abs(xvect(i)-xvect(i-1))/abs(xvect(i-1)-xvect(i-2)));
      p(i)=num/den;
      c(i)=abs(xvect(i+1)-xvect(i))/(abs(xvect(i)-xvect(i-1)))^p(i);
   end
end
