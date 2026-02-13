La richiesta del progetto riguardava la simulazione di uno sciame a partire da un singolo fotone, indagando la dipendenza dall'angolo di incidenza.

Lo sciame avviene a 20km e il rilevatore è posto a 4km dal suolo.
Per lo studio si utilizza un parametro chiamato cammino libero medio(lunghezza di radiazione) X0, questo parametro rappresenta anche l'unità di misura per la perdità di radiazione.
Infatti ad una profondità x l'energia residua sarà E(x)=exp[-x/X0], il valore di X= utilizzato in questo caso è X0=7\times10^4cm.
Le particelle coinvolte, elettroni e positroni sono cariche allora perdono energia per ionizzazione questo è il processo dominante, un fenomeno che a grande scala sembra continuo.
Mentre per i fotoni il processo dominante è la produzione di coppie elettrone-positrone con un libero cammino pari a 9X0/7.
Espresso il percorso interno in unità di lunghezze di radiazione t=x/X0 si ha che il numero di particelle sarà pari a n = 2^t.
Data l'energia delle perticelle incidenti E_0, l'energia di ogni particella sarà E(t)=E_0/2^t, avrò il massimo per E(t)=E_C.
Allora al momento della profondità massima dello sciame avremo t_max e quindi t_max=ln(E_0/E_C)/ln2. 

Per modellizzare questo fenomeno si utilizza il modello di Rossi.
Ogni interazione avviene dopo una lunghezza di interazione, ogni particella secondaria eredita metà energia della madre, per quanto riguarda le singole particelle:
l'interazione di un elettrone-positrone con energia E porta alla produzione di un fotone e un elettrone o positrone, mentre un fotone con energia E genera un elettrone ed un positrone.
Il processo si arresta quando l'energia degli elettroni scende sotto valore di energia critica E_c, l'energia quando la perdita per ionizzazione eguaglia quella per radiazione. 

L'utente può scegliere: l'energia iniziale tra 1 e 100 TeV, il passo di profondità "s" tra zero e uno (compreso), con un passo di 0.1 in 0.1.
Inoltre viene svolto uno studio statistico chiedendo all'utente il numero di eventi da simulare.

Il programma restituisce il profilo dei tre sciami a rispettivamente 0, 20 e 40 gradi. 
Lo studio statistico invece rappresenta un istogramma che ci da la densità rispetto al numero di hit per i tre eventi si è inoltre sviluppato un fit con uno studio lognormale, poiché ho un fenomeno di tipo stocastico, per verificare l'andamento. 
Il fit lognormale restituisce due parametri \mu e \sigma, il primo rappresenta la media della distribuzione del \ln(N) mentre il secondo la deviazione standard del \ln(N). Con il significato di media e fluttuazione dello sciame.

Un altro tipo di studio che è stato svolto è quello che mette in relazione il numero medio di hit al crescere dell'energia. Il tutto fittato con una retta, e il parametro \alpha della retta rappresentano quanto è lineare questa crescita.

L'ultimo grafico rappresenta come cambia il rapporto \sigma/\mu che misura quanto uno sciame oscilla da un evento all’altro, ovvero la stabilità dello sciame.
NB questo serve a mostrare come il modello di Rossi utilizzato sia limitante per il vero comportamento dello sciame.


