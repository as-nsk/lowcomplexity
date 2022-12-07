#ifndef DIVIZION2
#define DIVIZION2
#include "defines2.h"
#define NEW_SEQUENCE_ITEM 255
struct info{
	   byte method;//1 forward 2 reverse 3 forward func 4 rev func
	   s32bit from;
	   s32bit len;
	   s32bit seq_number;
	   };
int **tree;
struct list_store{
		 list *head;
		 list *tail;
		 };
list_store *low_level;//April 2003
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
			      //limit максимальный индекс начала вхождения
			      //size длина подпоследовательности
			      //buff подпоследовательность
			      //alf алфавит преобразования
			      //обе процедуры требуют вызова set_universe
			      virtual info find_forward_match(byte *buff,int size,int limit,byte *alf);
			      virtual info find_reverse_match(byte *buff,int size,int limit,byte *alf);
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
{info res={0,0,0,0};
 buff=buff;
 size=size;
 limit=limit;
 alf=alf;
 return res;
}
info abstract_divizion::find_reverse_match(byte *buff,int size,int limit,byte *alf)
{info res={0,0,0,0};
 buff=buff;
 size=size;
 limit=limit;
 alf=alf;
 return res;
}
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
 for(i=0;i<(NSIZE-1);i++)
    {
    citem=citem*ALPHABET_LEN+((int)buff[i]);
    if(tree[i+1][citem]==-1)tree[i+1][citem]=ent;
    }
 citem=citem*ALPHABET_LEN+((int)buff[i]);
 link_to_tail(citem,ent);
}
void divizion::add_vector(byte *buff,int ent,int max)
{int i,citem=0;
 if(max)
   {
   if(max<NSIZE)
     {
     for(i=0;i<max;i++)
	{
	citem=citem*ALPHABET_LEN+((int)buff[i]);
	if(tree[i+1][citem]==-1)tree[i+1][citem]=ent;
	 }
     }else add_vector(buff,ent);
   }
}


info divizion::find_forward_match(byte *buff,int size,int limit,byte *alf)
{info res={1,-1,0,0};
 if(alf[0]!=0)res.method=3;
 int citem=0;
 while((res.len<NSIZE)&&(tree[res.len][citem]!=-1)&&(res.len<size))//>=
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
	 curr=res.len;
	 from=hlp->from;
	 if(from<limit)
	   {
	   if((unisize-from)>size)maxlen=size;
          else maxlen=unisize-from;
	   while((alf[buff[curr]]==universe[from+curr])&&(curr<maxlen))curr++;
	   if(curr>max)
	     {
	     max=curr;
	     res.from=from;
	     }
	   }
	 hlp=hlp->next;
	 }
    res.len=max;
    }else res.len--;
   }else if(tree[res.len][citem]==-1)res.len--;
 return res;
}

info divizion::find_reverse_match(byte *buff,int size,int limit,byte *alf)
{info res={2,-1,0,0};
 if(alf[0]!=0)res.method=4;//method - inverse
 int citem=((int)alf[buff[0]]),mult=ALPHABET_LEN,i;
 i=1;
 while( (i<NSIZE)&&(size>=i) )//&&((tree[i][citem]+i)<limit)&&(size>=i) )
 {
// for(i=1;i<NSIZE;i++)//,i<limit
    if((tree[i][citem]!=-1)&&((tree[i][citem]+i)<=limit)&&(size>=i))
		{
        res.len=i;
        res.from=tree[i][citem];
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

//if(res.len>=NSIZE)
{
 if((low_level[citem].head!=NULL)&&(size>=NSIZE)&&(limit>=NSIZE))
   {int max=res.len,curr=0;
    list *hlp=low_level[citem].head;
    while(hlp!=NULL)
	 {int from;
	 curr=0;
	 from=hlp->from+NSIZE-1;
	 if(from<limit)
	   {
	   while((alf[buff[curr]]==universe[from-curr])&&(from>=curr)&&(curr<size))curr++;
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
{info res={2,0,1,0};
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
{info res={3,-1,1,0};
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
info clean_item={0,0,0,0};
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
 info result={0,0,0,0};
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
 info result={0,0,0,0};
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
