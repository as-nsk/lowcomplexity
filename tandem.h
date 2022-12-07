int findminrepunit(byte *seq,int seqsize)//find minimal rpeat unit in (tandem) repeated sequence
{
int j,minrepsize,i,flag=1;
minrepsize=seqsize; //no repeated units here
i=seqsize/2;

while(i>0)
{
j=0;
flag=1;
while(flag&&(j<seqsize-i))
	{
	if(seq[j]!=seq[j+i]) flag=0;//no coincindence,no repeat
	j++;
	}
if(flag) minrepsize=i;//fix minimal repeated size
i--;
}

return (minrepsize); //size of repeated unit
}

int revtandem=0,revpalindrome=0,dirtandem=0,dirclose=0;
char *legend2="atgc";
int	LongUniq=0,LongSSR=0;

int show_tandem_report(FILE *tandemfile, byte *buff0, info h0, int j, int limit, int hows)
	   {
	int i2,j2;

//	if(ALPHABET_LEN!=4) {printf("sorry,show_tandem_report is only for ATGC");exit(1);}

	if((limit<=h0.from-1) && ((h0.method==2)||(h0.method==4)) ) //find reversed repeat
		{
        fprintf(tandemfile,"Palindrome length=%i, pos=%i (lim=%i,from=%i,diff=%i), met=%i\n",h0.len,j, limit+hows, h0.from,(j-h0.from),h0.method);
		for(j2=0;j2<(limit-h0.from)+2*h0.len;j2++)
		  fprintf(tandemfile,"%c",legend2[buff0[h0.from-h0.len+j2]]);
          fprintf(tandemfile,"\n");
		for(j2=0;j2<h0.len;j2++)
		  fprintf(tandemfile,"<");
		for(j2=0;j2<(limit-h0.from);j2++)
		  fprintf(tandemfile,"-");
		for(j2=0;j2<h0.len;j2++)
		  fprintf(tandemfile,">");
          fprintf(tandemfile,"\n");

//		  i2=findminrepunit(buff0+hows+j,h0.len);
		  i2=findminrepunit(buff0+hows+(h0.from-h0.len),h0.len+h0.len-(j-h0.from));
          if(i2<h0.len)
				{
			    fprintf(tandemfile,"Reverse->Tandem - minimal rep.unit (%i bp *%.1f):",i2,(double)(h0.len)/(double)(i2));
			    revtandem++;
				}
		  else 
				{
				fprintf(tandemfile,"Reverse palindrome, no rep.unit (%i bp=%i): (",i2,(int)h0.len);
				revpalindrome++;
				}
		  for(j2=0;j2<i2;j2++)
	          fprintf(tandemfile,"%c",legend2[buff0[hows+j+j2]]);
          fprintf(tandemfile,")\n");
fprintf(tandemfile,"--- subroutine check -----------------\n");
		}

		if((h0.len>(limit-h0.from))&&(h0.len>1)&& ((h0.method==1)||(h0.method==3))) //find forward (tandem) repeat
{          
//show_item(h0,buff0+hows+j,j+hows);//(h0,buff0+hows+j,j+hows)
          fprintf(tandemfile,"Overlap length=%i, pos=%i (lim=%i,from=%i,diff=%i), met=%i\n",h0.len,j, limit, h0.from,(j-h0.from),h0.method);
//fprintf(tandemfile,"---------------\n");
	if((hows+j-(limit-h0.from+5))>0)
			for(j2=0;j2<(limit-h0.from+5);j2++)
			  fprintf(tandemfile,"%c",legend2[buff0[hows+j-(limit-h0.from+5)+j2]]);
          //fprintf(tandemfile,"[");
		  for(j2=0;j2<h0.len;j2++)
			  fprintf(tandemfile,"%c",legend2[buff0[hows+j+j2]]);
          //fprintf(tandemfile,"]");
		  for(j2=0;j2<5;j2++)
			  fprintf(tandemfile,"%c",legend2[buff0[hows+j+h0.len+j2]]);
          fprintf(tandemfile,"\n");
		  //mark coincindence
	if((hows+j-(limit-h0.from+5))>0)
			  for(j2=0;j2<(limit-h0.from+5);j2++) fprintf(tandemfile,"|");
		  for(j2=0;j2<h0.len;j2++)  fprintf(tandemfile,">");
		  for(j2=0;j2<5;j2++)  fprintf(tandemfile,"|");
          fprintf(tandemfile,"\n");
	if((hows+j-(limit-h0.from+5))>0)
		  for(j2=0;j2<5;j2++)
			  fprintf(tandemfile,"%c",legend2[buff0[h0.from-5+j2]]);
		  for(j2=0;j2<h0.len;j2++)
	          fprintf(tandemfile,"%c",legend2[buff0[h0.from+j2]]);
	if((hows+j-(limit-h0.from+5))>0)
		  for(j2=0;j2<5;j2++)
	          fprintf(tandemfile,"%c",legend2[buff0[h0.from+h0.len+j2]]);
          fprintf(tandemfile,"\n");
	if((hows+j-(limit-h0.from+5))>0)
		  for(j2=0;j2<5;j2++) fprintf(tandemfile,"-");
		  for(j2=0;j2<h0.len;j2++)  fprintf(tandemfile,"*");
		  for(j2=0;j2<5;j2++)  fprintf(tandemfile,"-");
          fprintf(tandemfile,"\n");

//		i2=findminrepunit(buff0+h0.from,h0.len);
		  i2=findminrepunit(buff0+h0.from,h0.len+(j-h0.from));

		if(i2<h0.len)
			{
			fprintf(tandemfile,"Tandem - minimal rep.unit (%i bp *%.1f):",i2,(double)(h0.len)/(double)(i2));
			dirtandem++;
			}
		else 
			{
			fprintf(tandemfile,"Artefact(Long tandem), rep.unit (%i bp=%i): (",i2,(int)h0.len);
			dirclose++;
			}

		for(j2=0;j2<i2;j2++)
	          fprintf(tandemfile,"%c",legend2[buff0[h0.from+j2]]);
          fprintf(tandemfile,")\n");
fprintf(tandemfile,"------------------------ subroutine check --\n");
		}

		return(0);
	   }

int show_long_repeats(FILE *tandemfile, byte *buff0, info h0, int j, int limit, int hows)
{
	int j2,i2;
          fprintf(tandemfile,"length=%i, pos=%i (lim=%i,from=%i,diff=%i), met=%i\n",h0.len,j, limit, h0.from,(j-h0.from),h0.method);
//		  show_item(h0,buff0+hows+j,j);//(h0,buff0+hows+j,j+hows)
		  for(j2=0;j2<h0.len;j2++)
			  fprintf(tandemfile,"%c",legend2[buff0[hows+j+j2]]);
          fprintf(tandemfile,"\n");
		 
		 if((h0.len>(limit-h0.from))&&(h0.len>1))
		 {
		 for(j2=0;j2<(h0.len+limit-h0.from);j2++)
			  fprintf(tandemfile,"%c",legend2[buff0[hows+h0.from+j2]]);
		 fprintf(tandemfile,"\n");
         fprintf(tandemfile,"Overlap length=%i, pos=%i (lim=%i,from=%i,diff=%i), met=%i\n",h0.len,j, limit, h0.from,(j-h0.from),h0.method);
//fprintf(tandemfile,"---------------\n");
	if((hows+j-(limit-h0.from+5))>0)
			for(j2=0;j2<(limit-h0.from+5);j2++)
			  fprintf(tandemfile,"%c",legend2[buff0[hows+j-(limit-h0.from+5)+j2]]);
          //fprintf(tandemfile,"[");
		  for(j2=0;j2<h0.len;j2++)
			  fprintf(tandemfile,"%c",legend2[buff0[hows+j+j2]]);
          //fprintf(tandemfile,"]");
		  for(j2=0;j2<5;j2++)
			  fprintf(tandemfile,"%c",legend2[buff0[hows+j+h0.len+j2]]);
          fprintf(tandemfile,"\n");
		  //mark coincindence
	if((hows+j-(limit-h0.from+5))>0)
			  for(j2=0;j2<(limit-h0.from+5);j2++) fprintf(tandemfile,"|");
		  for(j2=0;j2<h0.len;j2++)  fprintf(tandemfile,">");
		  for(j2=0;j2<5;j2++)  fprintf(tandemfile,"|");
          fprintf(tandemfile,"\n");
	if((hows+j-(limit-h0.from+5))>0)
		  for(j2=0;j2<5;j2++)
			  fprintf(tandemfile,"%c",legend2[buff0[h0.from-5+j2]]);
		  for(j2=0;j2<h0.len;j2++)
	          fprintf(tandemfile,"%c",legend2[buff0[h0.from+j2]]);
	if((hows+j-(limit-h0.from+5))>0)
		  for(j2=0;j2<5;j2++)
	          fprintf(tandemfile,"%c",legend2[buff0[h0.from+h0.len+j2]]);
          fprintf(tandemfile,"\n");
	if((hows+j-(limit-h0.from+5))>0)
		  for(j2=0;j2<5;j2++) fprintf(tandemfile,"-");
		  for(j2=0;j2<h0.len;j2++)  fprintf(tandemfile,"*");
		  for(j2=0;j2<5;j2++)  fprintf(tandemfile,"-");
          fprintf(tandemfile,"\n");
		 }

		i2=findminrepunit(buff0+hows+j,h0.len);
	j2=0;//-------additional check for hidden repeats
	while((j2<h0.len-1) && (i2>=h0.len))
		{
		i2=findminrepunit(buff0+hows+j,h0.len+j2);
		j2++;
		}

          if(i2<h0.len)
				{
				fprintf(tandemfile,"Inside long fragment->Tandem - minimal rep.unit (%i bp *%.1f) (+added %i):",i2,(double)(h0.len)/(double)(i2),j2);
				LongSSR++;
				}
		  else 
				{
				fprintf(tandemfile,"Long fragment - no simple repeat (%i): ",(int)h0.len);
				LongUniq++;
				}

		  if(i2<(int)h0.len)i2=(int)h0.len;
		  for(j2=0;j2<i2;j2++)        //(int)h0.len
	          fprintf(tandemfile,"%c",legend2[buff0[hows+j+j2]]);
          fprintf(tandemfile,"\n");

fprintf(tandemfile,"--------------------- subroutine check (long repeats) --\n");
return(0);
}