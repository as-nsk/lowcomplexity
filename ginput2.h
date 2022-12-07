#ifndef GINPUT2
#define GINPUT2
#include "defines2.h"
#define BUFF_SIZE 1024
#define END_OF_DATA -1
#define NEW_SEQUENCE 1000

char Seq_names[1000][20];//add April 2003
int Curr_seq_num=0;


class geninput {
		      char *buffer;
		      int *indexes;
		      int ALPHABET_LEN;
		      int CONTEXT_LEN;
		      u32bit **context;
		      void set_indexes(char *str);
	      protected:
			u32bit buffid,buffsize;
			void count_tail(u32bit tail,int len);
			virtual u32bit read_block(char *to,u32bit how);
	      public:
		     char legend[MAX_ALPHABET_LEN];
		     void input_init(char *s);
		     geninput();
		     int get_code();
		     int get_value(char i);
		     u32bit get_context();
		     int success_init();//return 1 if OK
		     int fill_context(int cont_len,u32bit **cont);
		     int get_alphabet_len();
		     ~geninput();
	      };
//-------------------------------------
geninput::geninput()
{
}
void geninput::input_init(char *s)
{
 ALPHABET_LEN=0;
 CONTEXT_LEN=0;
 indexes=new int[256];
 buffer=new char[BUFF_SIZE];
 if((indexes!=NULL)&&(buffer!=NULL))
   {
   buffid=buffsize=0;
   set_indexes(s);
   }
}
int geninput::get_alphabet_len()
{
 return ALPHABET_LEN;
}
int geninput::success_init()
{
 if((indexes!=NULL)&&(buffer!=NULL))return 1;
   else return 0;
}
int geninput::get_value(char i)
{
 return indexes[i];
}
u32bit geninput::read_block(char *to,u32bit how)
{
 to=to;
 how=how;
 return 0;
}
int geninput::get_code()
{int result=-1;
 int rem_exist=0;
 int item;
 do{
    while((result<0)&&(buffid<buffsize))
	 {
	 item=buffer[buffid++];
	 switch(item)
	       {
	       case '>':
	       case '<':rem_exist=2;
	       case '/':if(!rem_exist)rem_exist=1;//1-îáûêíîâåííûé êîììåíòàðèé
				   // 2-íîâàÿ ïîñëåäîâàòåëüíîñò
			result=-1;
			break;
	       case '\n':if(rem_exist==2)result=NEW_SEQUENCE;
			   else result=-1;
			 rem_exist=0;
			 break;
	       default:if(!rem_exist)result=indexes[item];
			 else result=-1;
	       }
	 }
    if(buffsize<=buffid)
      {
      buffsize=read_block(buffer,BUFF_SIZE);
      buffid=0;
      }else return result;
    }while(buffsize>buffid);
 if(result==-1)return END_OF_DATA;
   else return result;
}
/*
int geninput::get_code()
{int result=-1;
 int rem_exist=0;
 int item;
 int i;//,j;
 do{
    while((result<0)&&(buffid<buffsize))
	 {
	 item=buffer[buffid++];
	 switch(item)
	       {
	       case '>':
	       case '<':rem_exist=2;
	       case '/':
			   if(!rem_exist)
			rem_exist=1;//1-®¡ëª­®¢¥­­ë© ª®¬¬¥­â à¨©
				   // 2-­®¢ ï ¯®á«¥¤®¢ â¥«ì­®áâ
//-------------- add April 2003
//			if(Curr_seq_num<1000) //temporary limit - to change
//			{
//			i=0;
//			while((buffid<buffsize)&&(buffer[buffid]!='\n')&&(i<20))
//				{
//				Seq_names[Curr_seq_num][i]=buffer[buffid];
//				i++;buffid++;
//				}
//			  printf("normal %s \n",Seq_names[Curr_seq_num]);
//			if((buffsize<=buffid)&&(buffer[buffid]!='\n')&&(i<20))//read next if it is not enough
//		      {
//				if(rem_exist==2)result=NEW_SEQUENCE;
//				rem_exist=3;
//			}
//		      buffsize=read_block(buffer,BUFF_SIZE);
//		      buffid=0;
//			  printf("second read! %s \n",Seq_names[Curr_seq_num]);
//			  }//continue construct name (length 20)
//			while((buffid<buffsize)&&(buffer[buffid]!='\n')&&(i<20))
//				{
//				Seq_names[Curr_seq_num][i]=buffer[buffid];
//				i++;buffid++;
//			    if(i>=20)printf("again %s \n",Seq_names[Curr_seq_num]);
//				}
//			Seq_names[Curr_seq_num][i]='\n';
//			Curr_seq_num++;
			
//		if(buffer[buffid]=='\n') {rem_exist=0;result=NEW_SEQUENCE;printf("OK end string\n");}
			//for(j=i;j<20;j++) 
			//   Seq_names[Curr_seq_num-1][j]='-';
//			}
//------------------------			
			result=-1;
			break;
	       case '\n':if(rem_exist==2)result=NEW_SEQUENCE;
			   else result=-1;
			 rem_exist=0;
			 break;
	       default:if(!rem_exist)result=indexes[item];
			 else result=-1;
	       }
	 }
    if(buffsize<=buffid)
      {
      buffsize=read_block(buffer,BUFF_SIZE);
      buffid=0;
      }else return result;
    }while(buffsize>buffid);
 if(result==-1)return END_OF_DATA;
   else return result;
}
*/
void geninput::count_tail(u32bit tail,int len)
{int i;
 u32bit divizor=1;
 for(i=0;i<len;i++)divizor*=ALPHABET_LEN;
 tail%=divizor;
 for(i=len;i>0;i--)
    {
    divizor/=ALPHABET_LEN;
    context[i][tail]++;
    tail%=divizor;
    }
}
u32bit geninput::get_context()
{u32bit result=0,r1;
 for(int i=0;i<=CONTEXT_LEN;i++)
    {
    r1=get_code();
    if(r1==NEW_SEQUENCE)
      {
      //count_tail(result,i);
      result=0;
      i--;
      continue;
      }
    result=result*ALPHABET_LEN+r1;
    }
 return result;
}
void geninput::set_indexes(char *str)
{ALPHABET_LEN=0;
 int next_item_flag=1;
 for(int i=0;i<256;i++)indexes[i]=-1;
 while(*str)
      {
      switch(*str)
	    {
	    case '[':next_item_flag=0;
		     break;
	    case ']':if(!next_item_flag)
		       {
		       next_item_flag=1;
		       ALPHABET_LEN++;
		       }
		     break;
	    default:
		    if(indexes[*str]==-1)
		      {
		      indexes[lowcase(*str)]=ALPHABET_LEN;
		      indexes[upcase(*str)]=ALPHABET_LEN;
		      legend[ALPHABET_LEN]=*str;
		      if(next_item_flag)ALPHABET_LEN++;
		      }
	    }
      str++;
      }
 legend[ALPHABET_LEN]=0;
}
int geninput::fill_context(int cont_len,u32bit **cont)
{u32bit index=0;
 int i,currcode;
 context=cont;
 u32bit how_context=1;
 for(i=0;i<(cont_len+1);i++)how_context*=ALPHABET_LEN;
 CONTEXT_LEN=cont_len;
 buffer=new char[BUFF_SIZE];
 if(buffer)
   {
   buffsize=buffid=0;
   index=get_context();
   currcode=get_code();
   while(currcode!=END_OF_DATA)
	{
	context[CONTEXT_LEN+1][index]++;
	index=(index*ALPHABET_LEN+currcode)%how_context;
	currcode=get_code();
	if(currcode==NEW_SEQUENCE)
	  {
	  //count_tail(index,CONTEXT_LEN+1);
	  context[CONTEXT_LEN+1][index]++;
	  context[CONTEXT_LEN][index%(how_context/ALPHABET_LEN)]++;
	  index=get_context();
	  currcode=get_code();
	  }
	}
   if(index<how_context)
     {
     context[CONTEXT_LEN+1][index]++;
     context[CONTEXT_LEN][index%(how_context/ALPHABET_LEN)]++;
     //count_tail(index,CONTEXT_LEN+1);
     }
   delete buffer;
   buffer=NULL;
   return 1;
   }else return 0;
}
geninput::~geninput()
{
 if(indexes)delete indexes;
 if(buffer)delete buffer;
}
//-------------------------------------
class genfileinput:public geninput{
		FILE *fin;
		protected:
		virtual u32bit read_block(char *to,u32bit how);
		public:
		       genfileinput();
		       void set_file(FILE *f);
		};
//-------------------------------------
genfileinput::genfileinput()
{
}
u32bit genfileinput::read_block(char *to,u32bit how)
{
 return fread(to,sizeof(char),how,fin);
}
void genfileinput::set_file(FILE *f)
{
 fin=f;
 buffid=buffsize=0;
}
//-------------------------------------
#endif