#include "tempo.h"
#include <stdio.h>
#include <string.h>

 /* retorna o tempo em segundos desde 1/jan/1970*/
long time70(void)
{
  return((long) time(0));
}

/* imprime o tempo em segundos dividido em dias, horas, minutos, seg.*/
void printtime(long t)
{
  long dias,horas,minutos,segundos;
  if (t<1) {printf("<1 segundo");return;}
  minutos = (long) t/60;
  segundos = t - minutos*60;
  horas = (long) minutos/60;
  minutos = minutos-horas*60;
  dias = (long) horas/24;
  horas = horas - dias*24;
  if (dias==1) printf("%ld dia",dias);
  if (dias>1) printf("%ld dias",dias);
  if (dias && (horas||minutos||segundos)) printf(", ");
  if (horas==1) printf("%ld hora",horas);
  if (horas>1) printf("%ld horas",horas);
  if (horas && (minutos||segundos)) printf(", ");
  if (minutos==1) printf("%ld minuto",minutos);
  if (minutos>1) printf("%ld minutos",minutos);
  if (minutos && segundos) printf(", ");
  if (segundos==1) printf("%ld segundo",segundos);
  if (segundos>1) printf("%ld segundos",segundos);
}

/* imprime o tempo em segundos dividido em dias, horas, minutos, seg.*/
void sprinttime(char *s,long t)
{
  long dias,horas,minutos,segundos;
  char aux[100];
  if (t<1) {sprintf(s,"%s","<1 segundo");return;}
  minutos = (long) t/60;
  segundos = t - minutos*60;
  horas = (long) minutos/60;
  minutos = minutos-horas*60;
  dias = (long) horas/24;
  horas = horas - dias*24;
  s[0]='\0';
  if (dias==1) {sprintf(aux,"%ld dia",dias);strcat(s,aux);}
  if (dias>1) {sprintf(aux,"%ld dias",dias);strcat(s,aux);}
  if (dias && (horas||minutos||segundos)) strcat(s,", ");
  if (horas==1) {sprintf(aux,"%ld hora",horas);strcat(s,aux);}
  if (horas>1) {sprintf(aux,"%ld horas",horas);strcat(s,aux);}
  if (horas && (minutos||segundos)) strcat(s,", ");
  if (minutos==1) {sprintf(aux,"%ld minuto",minutos);strcat(s,aux);}
  if (minutos>1) {sprintf(aux,"%ld minutos",minutos);strcat(s,aux);}
  if (minutos && segundos) strcat(s,", ");
  if (segundos==1) {sprintf(aux,"%ld segundo",segundos);strcat(s,aux);}
  if (segundos>1) {sprintf(aux,"%ld segundos",segundos);strcat(s,aux);}
}
/* imprime o tempo no formato dd,hh:mm:ss 
   mas caso não precise imprimir os termos de maior ordem, não imprime */
void shortprinttime(long t)
{
  long dias,horas,minutos,segundos;
  minutos = (long) t/60;
  segundos = t - minutos*60;
  horas = (long) minutos/60;
  minutos = minutos-horas*60;
  dias = (long) horas/24;
  horas = horas - dias*24;
  if (dias>0) printf("%1ldd,%2ldh:%2ldm:%2lds",dias,horas,minutos,segundos);
  else if (horas>0) printf("%1ldh:%2ldm:%2lds",horas,minutos,segundos);
  else if (minutos>0) printf("%2ldm:%2lds",minutos,segundos);
  else printf("%2lds",segundos);
}


