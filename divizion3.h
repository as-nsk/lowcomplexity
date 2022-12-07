#ifndef DIVIZION2
#define DIVIZION2
#include "defines2.h"
#define NEW_SEQUENCE_ITEM 255
struct info{
	   byte method;//1 forward 2 reverse 3 forward func 4 rev func
	   s32bit from;
	   s32bit len;
	   s32bit seq_number;
	   int scorenonexact;
	   };
int **tree;
struct list_store{
		 list *head;
		 list *tail;
		 };
list_store *low_level;//April 2003
//int LIMIT_NONEXACT=0;

//--------------------------------------
class abstract_divizion{
		       protected:
				 int **tree;
				 list_store *low_level;
				 int NSIZE;
				 int ALPHABET_LEN;
				 int success_init_flag;
				 byte *universe;
				 int unisize;
//				 int BASEunisize;
//				 int BYunisize;
				 int LIMIT_NONEXACT;
				 int link_to_tail(int citem,int ent);
				 void clear_low_level(int how);
				 virtual void init(int Nsize,int Alphabet_len){};
		       public:
			      abstract_divizion(){};
			      int success_init();
			      virtual void clear(){};
			      virtual void add_vector(byte *buff,int ent);
			      virtual void add_vector(byte *buff,int ent,int max);
			      void set_table(byte *buff,int how);
			      void set_universe(byte *uni,int size);
			      void set_nonexact(int nonexact);
			      //limit максимальный индекс начала вхождения
			      //size длина подпоследовательности
			      //buff подпоследовательность
			      //alf алфавит преобразования
			      //обе процедуры требуют вызова set_universe
			      virtual info find_forward_match(byte *buff,int size,int limit,byte *alf);
			      virtual info find_reverse_match(byte *buff,int size,int limit,byte *alf);

			      virtual info find_forward_match2(byte *buff,int size,byte *alf);
			      virtual info find_reverse_match2(byte *buff,int size,byte *alf);

				  virtual info find_forward_match3(byte *buff,int size, byte *alf, int fromfix, int minreplength);//all repeats located after 'fromfix' by length >=minreplength
			      virtual info find_reverse_match3(byte *buff,int size, byte *alf, int fromfix, int minreplength);//all repeats located after 'fromfix' by length >=minreplength

				  virtual int count_low_level(int how2, double* entropy);//November2003-Feb.2004
                  virtual int linguistic(int level);
				  virtual int linguistic_big(byte *buff,int how2, int *counteq);
				  virtual int linguistic_bigK(byte *buff,int how2, int *counteq, int Klimit);//new - Feb.2004

			      virtual void done(){};
			      ~abstract_divizion(){};
		       };
//--------------------------------------
int abstract_divizion::success_init()
{
 return success_init_flag;
}
void abstract_divizion::add_vector(byte *buff,int ent)// just seal
{
 buff=buff;
 ent=ent;
}
void abstract_divizion::add_vector(byte *buff,int ent,int max)
{
 buff=buff;
 ent=ent;
 max=max;
}
info abstract_divizion::find_forward_match(byte *buff,int size,int limit,byte *alf)
{info res={0,0,0,0,0};//5-scorenonexact
 buff=buff;
 size=size;
 limit=limit;
 alf=alf;
 return res;
}
info abstract_divizion::find_reverse_match(byte *buff,int size,int limit,byte *alf)
{info res={0,0,0,0,0};//5-scorenonxact
 buff=buff;
 size=size;
 limit=limit;
 alf=alf;
 return res;
}
//--------------------------------
info abstract_divizion::find_forward_match2(byte *buff,int size, byte *alf)
{info res={0,0,0,0,0};//5-scorenonexact
 buff=buff;
 size=size;
 alf=alf;
 return res;
}
info abstract_divizion::find_reverse_match2(byte *buff,int size, byte *alf)
{info res={0,0,0,0,0};//5-scorenonxact
 buff=buff;
 size=size;
 alf=alf;
 return res;
}
info abstract_divizion::find_forward_match3(byte *buff,int size, byte *alf, int fromfix, int minreplength)//all repeats located after 'fromfix' by length >=minreplength
{info res={0,0,0,0,0};
 buff=buff;
 size=size;
 alf=alf;
 fromfix=fromfix;
 minreplength=minreplength;
 return res;
}

info abstract_divizion::find_reverse_match3(byte *buff,int size, byte *alf, int fromfix, int minreplength)//all repeats located after 'fromfix' by length >=minreplength
{info res={0,0,0,0,0};
 buff=buff;
 size=size;
 alf=alf;
 fromfix=fromfix;
 minreplength=minreplength;
 return res;
}


//--------------------------------
void abstract_divizion::clear_low_level(int how)
{int i;
 tree[0][0]=0;
 for(i=0;i<how;i++)
    {
    if(low_level[i].head!=NULL)
      {list *next=low_level[i].head,*next1;
      while(next!=NULL)
	   {
	   next1=next->next;
	    delete next;
	    next=next1;
	   }
      }
    low_level[i].head=NULL;
    low_level[i].tail=NULL;
    }
}

//-------- New complexities ---------------------------------------------------------
int abstract_divizion::linguistic(int level)
{
	level=level;
	return 0;
}
int abstract_divizion::count_low_level(int how2, double* entropy)//seal
{
	how2=how2;
	entropy=entropy;
return 0;
}
int abstract_divizion::linguistic_big(byte *buff,int how2, int *counteq)//seal
{
	buff=buff;
	how2=how2;
	counteq=counteq;
return 0;
}
int abstract_divizion::linguistic_bigK(byte *buff,int how2, int *counteq, int Klimit)//new - Feb.2004
{
	buff=buff;
	how2=how2;
	counteq=counteq;
	Klimit=Klimit;
return 0;
}


int abstract_divizion::link_to_tail(int citem,int ent)
{list *hlp=new list;
 if(hlp!=NULL)
   {
   hlp->from=ent;
   hlp->next=NULL;
   if(low_level[citem].head!=NULL)
     {
     low_level[citem].tail->next=hlp;
     low_level[citem].tail=hlp;
     }else low_level[citem].head=low_level[citem].tail=hlp;
   return 1;
   }else
   {//hlp==NULL: not enought memory
	   printf("not enought memory - abstract_divizion::link_to_tail citem=%d from=ent=%d\n",citem,ent);

   if(low_level[citem].head!=NULL)
     {
     low_level[citem].tail->next=low_level[citem].head;
     low_level[citem].tail=low_level[citem].head;
     low_level[citem].head=low_level[citem].head->next;
     low_level[citem].tail->next=NULL;
     low_level[citem].tail->from=ent;
     }
   }
 return 0;
}
void abstract_divizion::set_universe(byte *uni,int size)
{
 universe=uni;
 unisize=size;
}
//void abstract_divizion::set_universe2(byte *uni,int BASEsize,int BYsize)
//{
// universe=uni;
// BASEunisize=BASEsize;
// BYunisize=BASEsize;
//}
void abstract_divizion::set_nonexact(int nonexact)
{
LIMIT_NONEXACT=nonexact;
}

void abstract_divizion::set_table(byte *buff,int how)
{int i,citem=0,j;
 for(i=0;i<=(how-NSIZE);i++)add_vector(buff+i,i);
 while(i<how)
    {
    for(j=0;(i+j)<how;i++)
       {
       citem=citem*ALPHABET_LEN+((int)buff[i+j]);
       if(tree[j+1][citem]==-1)tree[j+1][citem]=i;
       }
    i++;
    }
}
//--------------------------------------
class divizion:public abstract_divizion
	     {
	     protected:
		     virtual void init(int Nsize,int Alphabet_len);
	     public:
		    divizion(int Nsize,int Alphabet_len);
		    divizion(){};
		    void clear();
		    void add_vector(byte *buff,int ent);
		    void add_vector(byte *buff,int ent,int max);
		    info find_forward_match(byte *buff,int size,int limit,byte *alf);
		    info find_reverse_match(byte *buff,int size,int limit,byte *alf);

			info find_forward_match2(byte *buff,int size, byte *alf);
		    info find_reverse_match2(byte *buff,int size, byte *alf);
			info find_forward_match3(byte *buff,int size, byte *alf, int fromfix, int minreplength);//all repeats located after 'fromfix' by length >=minreplength
			info find_reverse_match3(byte *buff,int size, byte *alf, int fromfix, int minreplength);//all repeats located after 'fromfix' by length >=minreplength

			int count_low_level(int how2, double* entropy);//November2003-Feb.2004
            int linguistic(int level);
            int linguistic_big(byte *buff,int how2, int *counteq);
            int linguistic_bigK(byte *buff,int how2, int *counteq, int Klimit);//new - Feb.2004

			void done();
		    ~divizion(){};
	     };
//--------------------------------------
divizion::divizion(int Nsize,int Alphabet_len)
{
 init(Nsize,Alphabet_len);
}
void divizion::init(int Nsize,int Alphabet_len)
{int i,how=1;
 NSIZE=Nsize;
 ALPHABET_LEN=Alphabet_len;
 success_init_flag=1;//all ok
 tree=new int*[NSIZE];
 if(tree!=NULL)
   {
   for(i=0;i<NSIZE;i++)
      {
      tree[i]=new int[how];
      if(tree[i]==NULL)success_init_flag=0;
      how*=ALPHABET_LEN;
      }
   low_level=new list_store[how];
   if(low_level!=NULL)
     {
     for(i=0;i<how;i++)low_level[i].head=low_level[i].tail=NULL;
     clear();
     }else success_init_flag=0;
   }else success_init_flag=0;
}
void divizion::clear()
{int i,j,how=1;
 for(i=0;i<NSIZE;i++)
    {
    for(j=0;j<how;j++)tree[i][j]=-1;
    how*=ALPHABET_LEN;
    }
 clear_low_level(how);
}
void divizion::add_vector(byte *buff,int ent)
{int citem=0,i;
i=0;//check DNA15 September
while((i<NSIZE-1)&&((int)buff[i]<ALPHABET_LEN))// for(i=0;i<(NSIZE-1);i++)
    {
    citem=citem*ALPHABET_LEN+((int)buff[i]);
    if(tree[i+1][citem]==-1)tree[i+1][citem]=ent;
	i++;
    }
 if((i==NSIZE-1) && (int)buff[i]<ALPHABET_LEN )//check DNA15 September - verified October2003
 {
 citem=citem*ALPHABET_LEN+((int)buff[i]);
 link_to_tail(citem,ent);
 }
}
void divizion::add_vector(byte *buff,int ent,int max)
{int i,citem=0;
if((int)buff[0]<ALPHABET_LEN)//verified October
 if(max)// && ((int)buff[0]<ALPHABET_LEN))//check DNA15  September
   {
   if(max<NSIZE)
     {
    // for(i=0;i<max;i++)//check DNA15 September
	i=0;
	while((i<max)&&((int)buff[i]<ALPHABET_LEN))
	{
	citem=citem*ALPHABET_LEN+((int)buff[i]);
	if(tree[i+1][citem]==-1)tree[i+1][citem]=ent;
	i++;
	 }

     }
   else add_vector(buff,ent);
   }
}
info divizion::find_forward_match(byte *buff,int size,int limit,byte *alf)
{info res={1,-1,0,0,0};//5-scorenonexact
 if(alf[0]!=0)res.method=3;
 int citem=0;
 //check DNA15 ALPHABET_LEN
 //imperfect
	int score=0;
	int oldscore=0;//add September nonexact

 while((res.len<NSIZE)&&(tree[res.len][citem]!=-1)&&(res.len<size)&&((int)alf[buff[res.len]]<ALPHABET_LEN))//>=
      {
      res.from=tree[res.len][citem];
      citem=citem*ALPHABET_LEN+((int)alf[buff[res.len]]);
      res.len++;
      }
 //if( (int)alf[buff[res.len]]<ALPHABET_LEN ) res.len--;//October 2003 additional check for NNNN
 if(res.len>=NSIZE)
   {
   if(low_level[citem].head!=NULL)
   {int max=0,curr=0,maxlen;
    list *hlp=low_level[citem].head;
    while(hlp!=NULL)
	 {int from;
	 curr=res.len;
	 oldscore=0;
	 from=hlp->from;
	 if(from<limit)
	   {
	   if((unisize-from)>size)maxlen=size;
          else maxlen=unisize-from;
//------ imperfect
	  score=0;
	  while((score<=LIMIT_NONEXACT)&&(curr<maxlen))
	  {
	  if((alf[buff[curr]]!=universe[from+curr])||((int)alf[buff[curr]]>=ALPHABET_LEN)) score++;
	  curr++;
	  }
	  if(score>LIMIT_NONEXACT)//||(curr>maxlen)) //set back
	  {curr--;score--;}
      if(score>0)
	  {
	  while((curr>0)&&((alf[buff[curr]]!=universe[from+curr])||((int)alf[buff[curr]]>=ALPHABET_LEN)) ) 
		  curr--;
	  }
//check right non-exact flank - atgcNNN
	  //------ imperfect
	   if((curr-score)>(max-oldscore))
	     {
	     max=curr;
		 oldscore=score;
	     res.from=from;
	     }
	   }
	 hlp=hlp->next;
	 }
    res.len=max;
	res.scorenonexact=oldscore;
    }else res.len--;
   }else if(tree[res.len][citem]==-1)res.len--;
 return res;
}

info divizion::find_forward_match2(byte *buff,int size, byte *alf)
{info res={1,-1,0,0,0};//5-scorenonexact
 if(alf[0]!=0)res.method=3;
 int citem=0; //check DNA15 ALPHABET_LEN
	int score, oldscore;//add September nonexact

 while((res.len<NSIZE)&&(tree[res.len][citem]!=-1)&&(res.len<size)&&((int)alf[buff[res.len]]<ALPHABET_LEN))//>=
      {
      res.from=tree[res.len][citem];
      citem=citem*ALPHABET_LEN+((int)alf[buff[res.len]]);
      res.len++;
      }
 if(res.len>=NSIZE)
   {
   if(low_level[citem].head!=NULL)
   {int max=0,curr=0,maxlen;
    list *hlp=low_level[citem].head;
    while(hlp!=NULL)
	 {int from;
	 curr=res.len-1;
	 oldscore=0;
	 from=hlp->from;
//	 if(from<limit)
//	   {
	   if((unisize-from)>size)maxlen=size;
          else maxlen=unisize-from;
//------ imperfect
	  score=0;
	  while((score<=LIMIT_NONEXACT)&&(curr<maxlen))
	  {
	  if((alf[buff[curr]]!=universe[from+curr])||((int)alf[buff[curr]]>=ALPHABET_LEN)) score++;
	  curr++;
	  }
	  if(score>LIMIT_NONEXACT)//||(curr>maxlen)) //set back
	  {curr--;score--;}
      if(score>0)
	  {
	  while((curr>0)&&((alf[buff[curr]]!=universe[from+curr])||((int)alf[buff[curr]]>=ALPHABET_LEN)) ) 
		  curr--;
	  }
//check right non-exact flank - atgcNNN
	  //------ imperfect
	   if((curr-score)>(max-oldscore))
	     {
	     max=curr;
		 oldscore=score;
	     res.from=from;
	     }
//	   }
	 hlp=hlp->next;
	 }
    res.len=max;
	res.scorenonexact=oldscore;
    }else res.len--;	//---?
   }else if(tree[res.len][citem]==-1) res.len--;   //--------?
    // if( (int)alf[buff[res.len]]<ALPHABET_LEN ) res.len--;}//October 2003 additional check for NNNN
 return res;
}

//tomorrow
info divizion::find_forward_match3(byte *buff,int size, byte *alf, int fromfix, int minreplength)//all repeats located after 'fromfix' by length >=minreplength
{info res={1,-1,0,0,0};//5-scorenonexact
 if(alf[0]!=0)res.method=3;
 int citem=0; //check DNA15 ALPHABET_LEN
	int score;//, oldscore;//add September nonexact

 while((res.len<NSIZE)&&(tree[res.len][citem]!=-1)&&(res.len<size)&&((int)alf[buff[res.len]]<ALPHABET_LEN))//>=
      {
      res.from=tree[res.len][citem];
      citem=citem*ALPHABET_LEN+((int)alf[buff[res.len]]);
      res.len++;
      }

 if(res.len>=NSIZE)
   {
   
 if(low_level[citem].head!=NULL)
   {//int max=0,
	   int curr=0;
   int maxlen;
   int from;
   //positioning   
   from=-1;
   list *hlp=low_level[citem].head;
    while(hlp!=NULL && from<fromfix)
	 {
	 from=hlp->from;
	 if(from<fromfix) hlp=hlp->next;
	 }

	if(from>=fromfix)
    while(hlp!=NULL && res.len<minreplength)
	 {
	 curr=res.len-1;
	 from=hlp->from;
	 if((unisize-from)>size)maxlen=size;
          else maxlen=unisize-from;

	//------ imperfect
	  score=0;
	  while((score<=LIMIT_NONEXACT)&&(curr<maxlen))
	  {
	  if((alf[buff[curr]]!=universe[from+curr])||((int)alf[buff[curr]]>=ALPHABET_LEN)) score++;
	  curr++;
	  }
	  if(score>LIMIT_NONEXACT)//||(curr>maxlen)) //set back
	  {curr--;score--;}
      if(score>0)
	  {
	  while((curr>0)&&((alf[buff[curr]]!=universe[from+curr])||((int)alf[buff[curr]]>=ALPHABET_LEN)) ) 
		  curr--;
	  }
//check right non-exact flank - atgcNNN
	  //------ imperfect

	  if(curr>=minreplength)//stop after find one large repeat >=minreplength
	  {
	  res.from=from;
	  res.len=curr;
	  res.scorenonexact=score;
	  }
//	   }
	 hlp=hlp->next;
	 }
//    res.len=max;
//	res.scorenonexact=oldscore;
    }else res.len--;	//---?
 
 }
 else if(tree[res.len][citem]==-1) res.len--;   //--------?
    // if( (int)alf[buff[res.len]]<ALPHABET_LEN ) res.len--;}//October 2003 additional check for NNNN
 return res;
}


info divizion::find_reverse_match(byte *buff,int size,int limit,byte *alf)
{info res={2,-1,0,0,0};//5-scorenonexact
 if(alf[0]!=0)res.method=4;//method - inverse
 int citem=((int)alf[buff[0]]),mult=ALPHABET_LEN,i;
 
 if((int)alf[buff[0]]>=ALPHABET_LEN) return(res);
 i=1;//DNA 15
 while( (i<NSIZE)&&(size>=i)&&((int)alf[buff[i]]<ALPHABET_LEN) )//&&((int)alf[buff[res.len]]<ALPHABET_LEN) //&&((tree[i][citem]+i)<limit)&&(size>=i) )
 {
// for(i=1;i<NSIZE;i++)//,i<limit
    if((tree[i][citem]!=-1)&&((tree[i][citem]+i)<=limit)&&(size>=i))
		{
        res.len=i;
        res.from=tree[i][citem]+i-1;//+i-1
		}
	citem+=(((int)alf[buff[i]])*mult);
    mult*=ALPHABET_LEN;  
 i++;
 }
/*
 while((res.len<NSIZE)&&(tree[res.len][citem]!=-1)&&(res.len<size))//>=
      {
      res.from=tree[res.len][citem];
      citem=citem*ALPHABET_LEN+((int)alf[buff[res.len]]);
      res.len++;
      }
*/
//------crazy april 2003

if(res.len>=NSIZE-1)
{//DNA15
 if((low_level[citem].head!=NULL)&&(size>=NSIZE)&&(limit>=NSIZE)&&((int)alf[buff[res.len]]<ALPHABET_LEN))
   {int max=res.len,curr=0;
    list *hlp=low_level[citem].head;
    while(hlp!=NULL)
	 {int from;
	 curr=0;
	 from=hlp->from+NSIZE-1; //   ||||<-----||||
	 if(from<limit)
	   {
	   while((alf[buff[curr]]==universe[from-curr])&&(from>=curr)&&(curr<size)&&((int)alf[buff[curr]]<ALPHABET_LEN))curr++;
	   if(curr>max)
	     {
	     max=curr;
	     res.from=from;
	     }
	   }
	 hlp=hlp->next;
	 }
    res.len=max;
    }
}

 return res;
}

info divizion::find_reverse_match2(byte *buff,int size, byte *alf)//no 'limit'
{info res={2,-1,0,0,0};//5-scorenonexact
 if(alf[0]!=0)res.method=4;//method - inverse
 int citem=((int)alf[buff[0]]),mult=ALPHABET_LEN,i;
 
 if((int)alf[buff[0]]>=ALPHABET_LEN) return(res);
 i=1;//DNA 15
 while( (i<NSIZE)&&(size>=i)&&((int)alf[buff[i]]<ALPHABET_LEN) )//&&((int)alf[buff[res.len]]<ALPHABET_LEN) //&&((tree[i][citem]+i)<limit)&&(size>=i) )
 {
// for(i=1;i<NSIZE;i++)//,i<limit
    if((tree[i][citem]!=-1)&&(size>=i))	//&&((tree[i][citem]+i)<=limit)
		{
        res.len=i;
        res.from=tree[i][citem]+i-1;//+i-1
		}
	citem+=(((int)alf[buff[i]])*mult);
    mult*=ALPHABET_LEN;  
 i++;
 }

if(res.len>=NSIZE-1)
{//DNA15
 if((low_level[citem].head!=NULL)&&(size>=NSIZE)&&((int)alf[buff[res.len]]<ALPHABET_LEN))
   {int max=res.len,curr=0;
    list *hlp=low_level[citem].head;
    while(hlp!=NULL)
	 {int from;
	 curr=0;
	 from=hlp->from+NSIZE-1; //   ||||<-----||||
//	 if(from<limit)
//	   {
	   while((alf[buff[curr]]==universe[from-curr])&&(from>=curr)&&(curr<size)&&((int)alf[buff[curr]]<ALPHABET_LEN))curr++;
	   if(curr>max)
	     {
	     max=curr;
	     res.from=from;
	     }
//	   }
	 hlp=hlp->next;
	 }
    res.len=max;
    }
}

 return res;
}

info divizion::find_reverse_match3(byte *buff,int size, byte *alf, int fromfix, int minreplength)//all repeats>=minreplen, no 'limit'
{info res={2,-1,0,0,0};//5-scorenonexact
 if(alf[0]!=0)res.method=4;//method - inverse
 int citem=((int)alf[buff[0]]),mult=ALPHABET_LEN,i;
 
 if((int)alf[buff[0]]>=ALPHABET_LEN) return(res);
 i=1;//DNA 15
 while( (i<NSIZE)&&(size>=i)&&((int)alf[buff[i]]<ALPHABET_LEN) )//&&((int)alf[buff[res.len]]<ALPHABET_LEN) //&&((tree[i][citem]+i)<limit)&&(size>=i) )
 {
    if((tree[i][citem]!=-1)&&(size>=i))	//&&((tree[i][citem]+i)<=limit)
		{
        res.len=i;
        res.from=tree[i][citem]+i-1;//+i-1
		}
	citem+=(((int)alf[buff[i]])*mult);
    mult*=ALPHABET_LEN;  
 i++;
 }

if(res.len>=NSIZE-1)
{//DNA15
 if((low_level[citem].head!=NULL)&&(size>=NSIZE)&&((int)alf[buff[res.len]]<ALPHABET_LEN))
   {//int max=res.len,
   int curr=0;
    int from;
   from=-1;
   list *hlp=low_level[citem].head;
     while(hlp!=NULL && from+NSIZE-1<fromfix)
	 {
	 from=hlp->from;
	 hlp=hlp->next;
	 }
//printf("we are here from=%i\n",from);//debug
   if(from+NSIZE-1>=fromfix)
    while(hlp!=NULL && res.len<minreplength)
	 {
//		printf("!!!!!!!!we are here from=%i\n",from);//debug
	 curr=0;
	 from=hlp->from+NSIZE-1; //   ||||<-----||||
//	 if(from<limit)
//	   {
	   while((alf[buff[curr]]==universe[from-curr])&&(from>=curr)&&(curr<size)&&((int)alf[buff[curr]]<ALPHABET_LEN))curr++;
	   if(curr>=minreplength)
	     {
//	     max=curr;
		 res.len=curr;
	     res.from=from;
	     }
//	   }
	 hlp=hlp->next;
	 }
//    res.len=max;
    }
}

 return res;
}


int divizion::linguistic(int level)
{int i,j,how=1,level0=2,LC=0;
if(level>2) level0=level; else {printf("lingustic - bad level=%i\n",level);}
if(level0>NSIZE) level0=NSIZE;

 for(i=0;i<level0;i++)
    {
    for(j=0;j<how;j++)
		if(tree[i][j]>-1) LC++;
    how*=ALPHABET_LEN;
    }
return(LC);
}

int divizion::linguistic_big(byte *buff,int seqlen, int *counteq)//new - November 2003
{int i,j;
int count_word=0,LC=0;
int count1=0, count2=0,count3=0;
//int countinside=0;
int from1,from2;
double entrop=0;
int countlocal[10000];//10K-maximum for linguistic complexity
 //tree[0][0]=0;
 int how=1;
// for(i=0;i<NSIZE;i++)
 //   how*=ALPHABET_LEN;
 if(seqlen>10000) {printf("LC - bad length =%i >10K\n",seqlen);return(LC);}
 if(seqlen<NSIZE) {printf("LC - bad length =%i NSIZE=%i\n",seqlen,NSIZE);return(LC);}
 for(i=0;i<=seqlen;i++) {counteq[i]=0;countlocal[i]=0;}

 how=ALPHABET_LEN;
 for(i=1;i<NSIZE;i++)
    {
    for(j=0;j<how;j++)
		if(tree[i][j]>-1) {LC++;counteq[i]++;}
    how*=ALPHABET_LEN;
    }

int counthow=0;
int count_notword=0;
 for(i=0;i<how;i++)
    {
    if(low_level[i].head!=NULL)
      {
	  list *hlp=low_level[i].head,*hlp2;//,*next1;

      while(hlp!=NULL)
	   {
	 from1=hlp->from;
/*
	for(j=0;j<NSIZE;j++)
	 printf("%1i",(int)(buff[from1+j]) );
		printf("- from1=%i\n",from1);
*/
	 count1++;//first copy
//     count2=1;
     hlp2=hlp->next;
     while(hlp2!=NULL)
	 {
	 from2=hlp2->from;
//-------------------
//	 if(from2+NSIZE<seqlen)
	    //{
        j=0;
	    while((buff[from1+NSIZE+j]==buff[from2+NSIZE+j]) && (from2+NSIZE+j<seqlen) && (NSIZE+j<seqlen) )
			{
			if((int)(buff[from1+NSIZE+j])<ALPHABET_LEN) countlocal[NSIZE+j]++;
			j++;
			}
		count2++;//second etc. copy
		//}
     hlp2=hlp2->next;//inner cycle
	 }	

	 for(j=0;j<seqlen-NSIZE;j++)
		{
		if(countlocal[NSIZE+j]>0) 
			{counteq[NSIZE+1+j]++;count3++;}//+1 !!!!!!!!!!
		}
/*
	 for(j=0;j<=7;j++)//debug------------
		printf("%i ",countlocal[j]);
		printf("- countlocal\n");
	 for(j=0;j<=7;j++)//debug------------
		printf("%i ",counteq[j]);
		printf("- count3=%i from1=%i from2=%i\n",count3,from1,from2);
*/
	 for(j=0;j<=seqlen-NSIZE;j++)
		countlocal[NSIZE+j]=0;


       hlp=hlp->next;//big cycle
	   }

      count_word++;
      }
	else
		count_notword++;
	counthow++;
    }

 counteq[NSIZE]=count_word;
 LC+=count_word;
// 	 for(j=0;j<=seqlen-NSIZE;j++)
//		counteq[NSIZE+j]=0;
// if(count2<seqlen-NSIZE+1){printf("bug? count2=%i, seqlen=%i\n",count2,seqlen);}
// if(count2>seqlen-NSIZE+1){printf("bug? count2=%i, seqlen=%i\n",count2,seqlen);}

// printf("count1(%i)=%i count2=%i, count3(eq)=%i count_word=%i (how=%i)\n",NSIZE,count1,count2,count3,count_word,how);
// printf("count_notword=%i counthow=%i counteq[NSIZE]=%i\n",count_notword,counthow,counteq[NSIZE]);
return LC;
}

int divizion::linguistic_bigK(byte *buff,int seqlen, int *counteq, int Klimit)//new - Feb.2004
{int i,j;
int count_word=0,LC=0;
int count1=0, count2=0,count3=0;
int from1,from2;
double entrop=0;
int countlocal[100002];//10K-maximum for linguistic complexity
 int how=1;
 if(Klimit>100000) {printf("LCK - bad length =%i >10K\n",seqlen);return(LC);}
 if(seqlen<1) {printf("LCK - bad length =%i NSIZE=%i\n",seqlen,NSIZE);return(LC);}
 for(i=0;i<=Klimit;i++) {counteq[i]=0;countlocal[i]=0;}

 how=ALPHABET_LEN;
//for(i=1;i<NSIZE;i++)
//  {
//   for(j=0;j<how;j++)
//		if(tree[i][j]>-1) {LC++;counteq[i]++;}
//    how*=ALPHABET_LEN;
//  }
 if(seqlen<NSIZE) return(0);//return 0 for very short sequences
 if(Klimit<NSIZE)//not need count for all
	{
	 for(i=1;i<Klimit;i++)
		{
	    for(j=0;j<seqlen;j++)
			if(tree[i][j]>-1) {LC++;counteq[i]++;}
		how*=ALPHABET_LEN;
		}
	 return(LC);
	 }
//else //standard calculation
for(i=1;i<NSIZE;i++)
	{
    for(j=0;j<how;j++)
		if(tree[i][j]>-1) {LC++;counteq[i]++;}
    how*=ALPHABET_LEN;
	}

int counthow=0;
int count_notword=0;
 for(i=0;i<how;i++)
    {
    if(low_level[i].head!=NULL)
      {
	  list *hlp=low_level[i].head,*hlp2;//,*next1;

      while(hlp!=NULL)
	   {
	 from1=hlp->from;
	 count1++;//first copy

     hlp2=hlp->next;
     while(hlp2!=NULL)
	 {
	 from2=hlp2->from;
//-------------------
        j=0;
	    while((buff[from1+NSIZE+j]==buff[from2+NSIZE+j]) && (from2+NSIZE+j<seqlen) && (NSIZE+j<Klimit) )
			{
			if((int)(buff[from1+NSIZE+j])<ALPHABET_LEN) countlocal[NSIZE+j]++;
			j++;
			}
		count2++;//second etc. copy

	 hlp2=hlp2->next;//inner cycle
	 }	

	 for(j=0;j<Klimit-NSIZE;j++)
		{
		if(countlocal[NSIZE+j]>0) 
			{counteq[NSIZE+1+j]++;count3++;}//+1 !!!!!!!!!!
		}
	 for(j=0;j<Klimit-NSIZE;j++)
		countlocal[NSIZE+j]=0;


       hlp=hlp->next;//big cycle
	   }

      count_word++;
      }
	else
		count_notword++;
	counthow++;
    }

 counteq[NSIZE]=count_word;
 LC+=count_word;

return LC;
}

//--bad copy K
/*
int divizion::linguistic_bigK(byte *buff,int how2, int *counteq, int Klimit)//new - Feb.2004
{int i,j;
int count_word=0,LC=0;
int count2=0,count3=0;
//int countinside=0;
int from1,from2;
double entrop=0;
int countlocal[10000];//10K-maximum for linguistic complexity
 //tree[0][0]=0;
 int how=1;
// for(i=0;i<NSIZE;i++)
 //   how*=ALPHABET_LEN;
 if(Klimit>10000) {printf("LC - bad length =%i >10K\n",how2);return(LC);}
 if(how2<NSIZE) {printf("LC - bad length =%i NSIZE=%i\n",how2,NSIZE);return(LC);}
 for(i=0;i<Klimit;i++) {counteq[i]=0;countlocal[i]=0;}

 how=ALPHABET_LEN;
 
 if(Klimit<NSIZE)//not need count for all
	{
	 for(i=1;i<Klimit;i++)
		{
	    for(j=0;j<how;j++)
			if(tree[i][j]>-1) {LC++;counteq[i]++;}
		how*=ALPHABET_LEN;
		}
	 return(LC);
	 }

 for(i=1;i<NSIZE;i++)
    {
    for(j=0;j<how;j++)
		if(tree[i][j]>-1) {LC++;counteq[i]++;}
    how*=ALPHABET_LEN;
    }
 
 for(i=0;i<how;i++)
    {
    if(low_level[i].head!=NULL)
      {
	  list *hlp=low_level[i].head,*hlp2;//,*next1;
//	  countinside=1;
      while(hlp!=NULL)
	   {
//	   next1=next->next;
	 from1=hlp->from;

     hlp2=hlp->next;
     while(hlp2!=NULL)
	 {
	 from2=hlp2->from;
//-------------------
//    list *hlp=low_level[citem].head;
//    while(hlp!=NULL)
//	 int from;
//	 curr=res.len;
     
	 if(from2+NSIZE<how2)
	   {
       j=0;
	   while(buff[from1+NSIZE+j]==buff[from2+NSIZE+j] && (from2+NSIZE+j<how2) && (NSIZE+j<Klimit) )
	   {if(buff[from1+NSIZE+j]<ALPHABET_LEN) countlocal[NSIZE+j]++;j++;}
	//   while((alf[buff[curr]]==buff[from+curr])&&(curr<maxlen))curr++;
		}
//	 hlp=hlp->next;
//   delete next;
//	    next=next1;
	 count3++;
     hlp2=hlp2->next;//inner cycle
	 }	
		count2++;
//		countinside++;
	for(j=0;j<Klimit;j++)
		if(countlocal[NSIZE+j]>0) 
		{
		counteq[NSIZE+j]++;
		countlocal[NSIZE+j]=0;
		}

       hlp=hlp->next;//big cycle
	   }
//      if(countinside>0)
//	   entrop-= ( ((double)countinside )/(double)(how2-NSIZE) )*log( ((double)countinside)/( (double)(how2-NSIZE) ) );
      count_word++;
      }

    }

// if(count2<how2-NSIZE+1){printf("bug? count2=%i, how2=%i\n",count2,how2);}
 if(count2>how2-NSIZE+1){printf("bug? count2=%i, how2=%i\n",count2,how2);}
// if(count2<count_word){printf("bug? count2=%i, count_word=%i\n",count2,count_word);}
// if(count2<count_word){printf("bug? count2=%i, count_word=%i\n",count2,count_word);}
// if(count3>count2){printf("bug? count2=%i, count3=%i\n",count2,count3);}
 // entrop/=(double)(how-NSIZE);

return LC;
}
*/

int divizion::count_low_level(int how2, double* entropy)//new - November 2003
{int i;
int count_word=0;
int count2=0;
int countinside=0;
double entrop=0;
 //tree[0][0]=0;
 int how=1;
 for(i=0;i<NSIZE;i++)
    {
    //for(j=0;j<how;j++)tree[i][j]=-1;
    how*=ALPHABET_LEN;
    }

 for(i=0;i<how;i++)
    {
    if(low_level[i].head!=NULL)
      {
	  list *next=low_level[i].head,*next1;
	  countinside=1;
      while(next!=NULL)
	   {
	   next1=next->next;
//	    delete next;
	    next=next1;
		count2++;
		countinside++;
	   }
      if(countinside>0)
	   entrop-= ( ((double)countinside )/(double)(how2-NSIZE) )*log( ((double)countinside)/( (double)(how2-NSIZE) ) );
      count_word++;
      }
//    low_level[i].head=NULL;
//    low_level[i].tail=NULL;
    }

 if(count2<how2-NSIZE+1){printf("bug? count2=%i, how2=%i\n",count2,how2);}
 if(count2>how2-NSIZE+1){printf("bug? count2=%i, how2=%i\n",count2,how2);}
 if(count2<count_word){printf("bug? count2=%i, count_word=%i\n",count2,count_word);}
 // entrop/=(double)(how-NSIZE);
 *entropy=entrop;
return count_word;
}




void divizion::done()
{int i;
 clear();
 if(low_level!=NULL)delete low_level;
 for(i=0;i<NSIZE;i++)if(tree[i]!=NULL)delete tree[i];
 if(tree!=NULL)delete tree;
}



//--------------------------------------
class multidivizion:public abstract_divizion
	     {
	     private:
		     int total;
		     int *scaler;
		     int *ready_func;
		     int *try_func;
		     int *try1_func;
	     protected:
		    void init(int Nsize,int Alphabet_len);
		    int compare(byte *uni,byte* buff,int limit,int incr);
	     public:
		    multidivizion(int Nsize,int Alphabet_len);
		    multidivizion(){};
		    void add_vector(byte *buff,int ent);
		    void clear();
		    void add_vector(byte *buff,int ent,int max);
		    //limit максимальный индекс начала вхождения
		    //size длина подпоследовательности
		    //buff подпоследовательность
		    //alf алфавит преобразования
		    //обе процедуры требуют вызова set_universe
		    info find_forward_match(byte *buff,int size,int limit,byte *alf);
		    info find_reverse_match(byte *buff,int size,int limit,byte *alf);
		    void done();
		    ~multidivizion();
		     };
//--------------------------------------
multidivizion::multidivizion(int Nsize,int Alphabet_len)
{
 init(Nsize,Alphabet_len);
}
void multidivizion::init(int Nsize,int Alphabet_len)
{int i,multiplier;
 NSIZE=Nsize;
 ALPHABET_LEN=Alphabet_len;
 success_init_flag=1;//all ok
 tree=new int*[NSIZE-1];
 scaler=new int[NSIZE];//scaler[0] consist summ of all levels
 ready_func=new int[ALPHABET_LEN];
 try_func=new int[ALPHABET_LEN];
 try1_func=new int[ALPHABET_LEN];
 if((tree!=NULL)&&(scaler!=NULL))
   {
   scaler[0]=1;
   multiplier=2;
   for(i=1;i<NSIZE;i++)
      {
      scaler[i]=scaler[i-1]*multiplier;
      multiplier*=2;
      scaler[0]+=scaler[i];
      }
   tree[0]=new int;
   tree[0][0]=0;
   for(i=1;i<NSIZE-1;i++)
      {
      tree[i]=new int[scaler[i]];
      if(tree[i]==NULL)success_init_flag=0;
      }
   low_level=new list_store[scaler[NSIZE-1]];
   if(low_level!=NULL)
     {
     for(i=0;i<scaler[NSIZE-1];i++)low_level[i].head=low_level[i].tail=NULL;
     clear();
     }else success_init_flag=0;
   }else success_init_flag=0;
}
void multidivizion::clear()
{int i,j;
 tree[0][0]=0;
 for(i=1;i<NSIZE-1;i++)
    {
    for(j=0;j<scaler[i];j++)tree[i][j]=-1;
    }
 clear_low_level(scaler[NSIZE-1]);
}
void multidivizion::add_vector(byte *buff,int ent)
{
 add_vector(buff,ent,NSIZE);
}
void multidivizion::add_vector(byte *buff,int ent,int max)
{int i,j,maximal,citem=0;
 if(max>NSIZE)maximal=NSIZE;
   else maximal=max;
 for(i=0;i<maximal;i++)
    {
    for(j=0;j<i;j++)
       {
       citem*=2;
       if(buff[i]==buff[j])citem++;
       }
    if((i<(NSIZE-1))&&(tree[i][citem]==-1))tree[i][citem]=ent;
    }
 if(max>=NSIZE)link_to_tail(citem,ent);
}
int multidivizion::compare(byte *uni,byte* buff,int limit,int incr)
{int i,id;
 for(i=0;i<ALPHABET_LEN;i++)try_func[i]=try1_func[i]=255;
 for(i=0,id=0;i<limit;i++,id+=incr)
    {
    if(try1_func[uni[id]]==255)
      {
      if(try_func[buff[i]]==255)
	{
	try1_func[uni[id]]=buff[i];
	try_func[buff[i]]=uni[id];
	}else break;
      }else
       {
       if(try_func[buff[i]]==255)break;
	 else
	 {
	 if(try_func[buff[i]]!=uni[id])break;
	 }
       }
    }
 return i;
}
info multidivizion::find_reverse_match(byte *buff,int size,int limit,byte *alf)
{info res={2,0,1,0,0};//score non-exact
 int citem=0,mult=1,i,j,maxlen;
 for(i=0;i<(NSIZE-1);i++)
    {
    for(j=0;j<i;j++)
       {
       if(buff[i]==buff[j])citem+=mult;
       mult*=2;
       }
    if((tree[i][citem]!=-1)&&((tree[i][citem]+i+1)<=limit)&&(size>=(i+1)))
      {
      res.len=i+1;
      res.from=tree[i][citem];
      }
    }
 for(j=0;j<i;j++)
    {
    if(buff[i]==buff[j])citem+=mult;
    mult*=2;
    }
 if((low_level[citem].head!=NULL)&&(size>=NSIZE))
   {int max=0,curr=0;
    res.len++;
    list *hlp=low_level[citem].head;
    while(hlp!=NULL)
	 {int from;
	 curr=res.len;
	 from=hlp->from+curr-1;
	 if(from<limit)
	   {
	   if(size>from)maxlen=from;
	     else maxlen=size;
	   curr=compare(universe+from,buff,maxlen,-1);
	   if(curr>max)
	     {
	     for(i=0;i<ALPHABET_LEN;i++)ready_func[i]=try_func[i];
	     max=curr;
	     res.from=from-curr+1;
	     }
	   }
	 hlp=hlp->next;
	 }
    res.len=max;
    }
    else
    {
    res.len=compare(universe+res.from+res.len-1,buff,res.len,-1);
    for(i=0;i<ALPHABET_LEN;i++)ready_func[i]=try_func[i];
    }
 for(i=0;i<ALPHABET_LEN;i++)alf[i]=ready_func[i];
 return res;
}
info multidivizion::find_forward_match(byte *buff,int size,int limit,byte *alf)
{info res={3,-1,1,0,0};//5-scorenonexact
 int citem=0,i;
 while(res.len<NSIZE)
      {
      if(tree[res.len-1][citem]!=-1)
	{
	res.from=tree[res.len-1][citem];
	for(i=0;i<res.len;i++)
	   {
	   citem*=2;
	   if(buff[res.len]==buff[i])citem++;
	   }
      if(res.len>=size)break;
      res.len++;
      }else
      {
      res.len--;
      break;
      }
      }
 if((res.len>=NSIZE)&&(low_level[citem].head!=NULL))
   {int max=0,curr=0,maxlen;
    list *hlp=low_level[citem].head;
    while(hlp!=NULL)
	 {int from;
	 curr=res.len;
	 from=hlp->from;
	 if(from<limit)
	   {
	   if((unisize-from)<size)maxlen=unisize-from;
	     else maxlen=size;
	   curr=compare(universe+from,buff,maxlen,1);
	   if(curr>max)
	     {
	     max=curr;
	     res.from=from;
	     for(i=0;i<ALPHABET_LEN;i++)ready_func[i]=try_func[i];
	     }
	   }
	 hlp=hlp->next;
	 }
    res.len=max;
    }else
     {
     res.len=compare(universe+res.from,buff,res.len,1);
     for(i=0;i<ALPHABET_LEN;i++)ready_func[i]=try_func[i];
     }
 for(i=0;i<ALPHABET_LEN;i++)alf[i]=ready_func[i];
 return res;
}
void multidivizion::done()
{int i;
 clear();
 if(low_level!=NULL)delete low_level;
 for(i=0;i<NSIZE-1;i++)if(tree[i]!=NULL)delete tree[i];
 if(tree!=NULL)delete tree;
 delete scaler;
 delete ready_func;
 delete try_func;
 delete try1_func;
}
multidivizion::~multidivizion()
{
}
//--------------------------------------
void clean_list(list *it)
{list *hlp;
 while(it!=NULL)
	  {
	  hlp=it->next;
	  delete it;
	  it=hlp;
	  }
}
info clean_item={0,0,0,0,0};
struct list_info{
		info dt;
		list_info* next;
		};
struct grammer{
	      list *enter;
	      grammer *son;
	      };
class slide_divizion:public abstract_divizion
	{
	private:
		grammer f_head;
		grammer r_head;
		int ft_vector;
		grammer* add_new_sons();
		void link_entrance(grammer *place, int ent);
		void add_forward_vector(byte *buff,int ent,int how);
		void add_reverse_vector(byte *buff,int ent,int how);
		info get_flist_match(byte *buff,int size,list *lst,byte *alfa);
		info get_rlist_match(byte *buff,int size,list *lst,byte *alfa);
		list* get_below_list(list *it,int limit);
		void shift_ftree();
		void shift_rtree();
		void recurse_deleting_sons(grammer *it);
	protected:
		void init(int Nsize,int Alphabet_len);
	public:
		slide_divizion(int Nsize,int Alphabet_len);
		slide_divizion(){};
		virtual void add_vector(byte *buff,int ent,int max);
		info find_forward_match(byte *buff,int size,int limit,byte *alfa);
		info find_reverse_match(byte *buff,int size,int limit,byte *alfa);
		void slide_tree(int FRAME_SIZE);
		void clear();
		void done();
		~slide_divizion(){};
	};
//--------------------------------------------------------
slide_divizion::slide_divizion(int Nsize,int Alphabet_len)
{
 init(Nsize,Alphabet_len);
}
void slide_divizion::init(int Nsize,int Alphabet_len)
{
 NSIZE=Nsize;
 ALPHABET_LEN=Alphabet_len;
 ft_vector=0;
 f_head.enter=NULL;
 f_head.son=NULL;
 r_head=f_head;
 success_init_flag=1;
}
//-------------- essential functions ---------------------
grammer* slide_divizion::add_new_sons()
{grammer *result=new grammer;
 int i;
 result=new grammer[ALPHABET_LEN];
 for(i=0;i<ALPHABET_LEN;i++)
	{
	result[i].enter=NULL;
	result[i].son=NULL;
	}
 return result;
}
void slide_divizion::link_entrance(grammer *place, int ent)
{list *i=new list;
 if(i!=NULL)
	{
	if((place->enter==NULL)||(place->enter->from<ent))
		{
		i->from=ent;
		i->next=place->enter;
		place->enter=i;
		}else delete i;
	}

}
void slide_divizion::add_forward_vector(byte *buff,int ent,int how)
{int i,max,currsym;
 grammer *cp;
 max=(NSIZE<how)?NSIZE:how;
 if(how)
	{
	if(f_head.son==NULL)f_head.son=add_new_sons();
	cp=&f_head;
	for(i=0;i<max;i++)
		{
		currsym=buff[i];
		cp=cp->son+currsym;
		link_entrance(cp,ent);
		if(cp->son==NULL)cp->son=add_new_sons();
		}
	}
}
void slide_divizion::add_reverse_vector(byte *buff,int ent,int how)
{int i,max,currsym;
 grammer *cp;
 max=(NSIZE<(ent+1))?NSIZE:(ent+1);
 if(how)
	{
	if(r_head.son==NULL)r_head.son=add_new_sons();
	cp=&r_head;
	for(i=0;i<max;i++)
		{
		currsym=buff[-i];
		cp=cp->son+currsym;
		link_entrance(cp,ent);
		if(cp->son==NULL)cp->son=add_new_sons();
		}
	}
}
void slide_divizion::add_vector(byte *buff,int ent,int max)
{
 add_forward_vector(buff,ent,max);
 add_reverse_vector(buff,ent,max);
}
//------- end of essential functions ---------------------
//-------------- searching functions ---------------------
info slide_divizion::get_flist_match(byte *buff,int size,list *lst,byte *alfa)
{int i;
 int maxsize;
 info result={0,0,0,0,0};//5
 while(lst!=NULL)
	{
	i=0;
	maxsize=(size<(unisize-lst->from))?size:unisize-lst->from;
	while((alfa[buff[i]]==universe[lst->from+i])&&(i<maxsize))i++;
	if(result.len<i)
		{
		if(alfa[0]==0)result.method=1;
			else result.method=3;
		result.len=i;
		result.from=lst->from;
		}
	lst=lst->next;
	}
 return result;
}
info slide_divizion::get_rlist_match(byte *buff,int size,list *lst,byte *alfa)
{int i;
 int maxsize;
 info result={0,0,0,0,0};//5
 while(lst!=NULL)
	{
	i=0;
	maxsize=(size<(lst->from+1-ft_vector))?size:lst->from+1-ft_vector;
	while((alfa[buff[i]]==universe[lst->from-i])&&(i<maxsize))i++;
	if(result.len<i)
		{
		if(alfa[0]==0)result.method=2;
			else result.method=4;
		result.len=i;
		result.from=lst->from-i+1;
		}
	lst=lst->next;
	}
 return result;
}
list* slide_divizion::get_below_list(list *it,int limit)
{
	while((it!=NULL)&&(it->from>=limit))it=it->next;
	return it;
}
info slide_divizion::find_forward_match(byte *buff,int size,int limit,byte *alfa)
{int i=0;
 grammer *cp=&f_head;
 if(cp->son!=NULL)
	{
	while(get_below_list(cp->son[alfa[buff[i]]].enter,limit)!=NULL)
		{
		if((i>=size)||(i>=NSIZE))break;
		cp=cp->son+alfa[buff[i]];
		i++;
		}
	return get_flist_match(buff,size,get_below_list(cp->enter,limit),alfa);
	}else return clean_item;
}
info slide_divizion::find_reverse_match(byte *buff,int size,int limit,byte *alfa)
{int i=0;
 grammer *cp=&r_head;
 if(cp->son!=NULL)
	{
	while(get_below_list(cp->son[alfa[buff[i]]].enter,limit)!=NULL)
		{
		if((i>=size)||(i>=NSIZE))break;
		cp=cp->son+alfa[buff[i]];
		i++;
		}
	return get_rlist_match(buff,size,get_below_list(cp->enter,limit),alfa);
	}else return clean_item;
}
//------- end of searching functions ---------------------
void slide_divizion::shift_ftree()
{grammer *cp=&f_head;
 int i=0;
 list *work;
 while(cp->son[universe[ft_vector+i]].enter!=NULL)
	{	//tail cutting
	work=cp->son[universe[ft_vector+i]].enter;
	if(work!=NULL)
		{
		if(work->next!=NULL)
			{
			while(work->next->next!=NULL)work=work->next;
			delete work->next;//there is would be comparsing work->from==ft_vector
			work->next=NULL;
			}else
			{
			delete work;
			cp->son[universe[ft_vector+i]].enter=NULL;
			}
		}
	cp=cp->son+universe[ft_vector+i];
	i++;
	}
}
void slide_divizion::shift_rtree()
{grammer *cp=&r_head;
 int i=0,j;
 list *work;
 for(i=0;i<NSIZE;i++)
	{
	cp=&r_head;
	for(j=i;j>=0;j--)cp=cp->son+universe[ft_vector+j];
	work=cp->enter;//tail cutting
	if(work!=NULL)
	  {
	  if(work->next!=NULL)
	    {
	    while(work->next->next!=NULL)work=work->next;
	    delete work->next;
	    work->next=NULL;
	    }else
	    {
	    delete work;
	    cp->enter=NULL;
	    }
	}
 }
}
void slide_divizion::slide_tree(int FRAME_SIZE)
{int offset=ft_vector+FRAME_SIZE;
 int i;
 shift_ftree();
 for(i=0;i<NSIZE;i++)add_forward_vector(universe+offset-i,offset-i,i+1);
 shift_rtree();
 add_reverse_vector(universe+offset,offset,offset+1);
 ft_vector++;
}
void slide_divizion::recurse_deleting_sons(grammer *it)
{int i;
 clean_list(it->enter);
 if(it->son!=NULL)
   {
   for(i=0;i<ALPHABET_LEN;i++)
	  {
	  recurse_deleting_sons(it->son+i);
	  }
   delete it->son;
   }
}
void slide_divizion::clear()
{
 recurse_deleting_sons(&f_head);
 recurse_deleting_sons(&r_head);
 ft_vector=0;
 f_head.enter=NULL;
 f_head.son=NULL;
 r_head=f_head;
}
void slide_divizion::done()
{
 clear();
}
//--------------------------------------
#endif DIVIZION2
