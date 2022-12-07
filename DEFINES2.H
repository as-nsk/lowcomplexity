#ifndef DEFINES2
#define DEFINES2
#include <stdio.h>
#define MAX_CONTEXT_LEN 12
#define MIN_CONTEXT_LEN 2
#define DEFAULT_CONTEXT_LEN 3
#define MAX_ALPHABET_LEN 30
#define MIN_ALPHABET_LEN 2
#define PI 3.14159265358979
#define SWITCH_CHAR '-'
typedef long double my_real;
typedef int u32bit;//unsigned 32 bit size integer
typedef int s32bit;//signed 32 bit size integer
typedef unsigned char byte;
struct list{//список вхождений
	   s32bit from;
	   list *next;
	   };
list* add_to_list(list *tail,s32bit item)
{list *addition=new list;
 if(addition)
   {
   addition->next=NULL;
   addition->from=item;
   if(tail!=NULL)tail->next=addition;
   }
 return addition;
}
char upcase(char i)
{
 if((i>='a')&&(i<='z'))return (i&223);
   else return i;
}
char lowcase(char i)
{
 if((i>='A')&&(i<='Z'))return (i|32);
   else return i;
}
u32bit renumeric(u32bit it,int len,int ALPHABET_LEN)
{u32bit i;
 u32bit result=0;
 for(i=0;i<len;i++)
    {
    result=result*ALPHABET_LEN+(it%ALPHABET_LEN);
    it/=ALPHABET_LEN;
    }
 return result;
}
typedef unsigned char byte;
class exact{
		 int how;
		 byte *presence;
		 byte *item;
		 int first_time;
	   public:
		  exact(int h);
		  void reset();
		  byte* get_next();
		  ~exact();
	   };
exact::exact(int h)
{if((h<=0)||(h>=100))h=5;
 how=h;
 presence=new byte[how];
 item=new byte[how];
 reset();
}
void exact::reset()
{int i;
 for(i=0;i<how;i++)
    {
    item[i]=how-1-i;
    presence[i]=1;
    }
 first_time=1;
}
byte* exact::get_next()
{int digit,i,j=0;
 if(!first_time)
   {
 do
  {
  for(i=0;i<how;i++)presence[i]=0;
  i=1;
  digit=-1;
  while(i&&(digit<how-1))
       {
       digit++;
       item[digit]++;
       if(item[digit]>=how)item[digit]=0;
	  else i=0;
//       presence[item[digit]]++;
       }
  j=0;
  for(i=0;i<how;i++)
     {
     if(presence[item[i]])
       {
       j=1;
       break;
       }else presence[item[i]]++;
     }
  }while(j);
  }
 first_time=0;
 return item;
}
exact::~exact()
{
 delete item;
 delete presence;
}
unsigned long factorial(unsigned long i)
{
 if(i)return i*factorial(i-1);
   else return 1;
}
void make_original_name(char *base,char *mask,char *result)
{int i=0,j=0,k;
 while(mask[i])
      {
      switch(mask[i])
	    {
	    case '?':result[j]=base[j];
		     j++;
		     break;
	    case '*':k=0;
		     while((base[k]!='.')&&(base[k]))
			  {
			  result[j++]=base[k++];
			  }
		     break;
	    default:result[j]=mask[i];
		    j++;
	    }
      i++;
      }
 result[j]=0;
}
long filesize(FILE *stream)
{long curpos, length;
 curpos = ftell(stream);
 fseek(stream, 0L, SEEK_END);
 length = ftell(stream);
 fseek(stream, curpos, SEEK_SET);
 return length;
}


#endif
