#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ginput3.h"
#include "divizion3.h"
#include "factorial.h"
#define SEQUENCE_LEN 100000000 //10Mb
#define MAX_FRAME_LEN SEQUENCE_LEN
#define MIN_FRAME_LEN 1
#define DEFAULT_FRAME_LEN 50
#define MIN_SLIDE_STEP 1
#define MAX_N_SPEC_STAT 100 //100
#define MAX_N_SPEC_STAT2 500//2000
#define MaxNumberofSeqs 45000 //1024 - 43000 for splice sites
#define MAX_SEQUENCIES 45000
#define MAX_METHODS 5//6 - reserved for LCbig
#define HOW_KEYS 8

//char* task_path;

// double ln1_2=1.44269504088896; //1/ln2 to calculate log(bybase2) via standard ln
//char default_alphabet[30]="atgc";
//char protein_alphabet[21]="ACDEFGHIKLMNPQRSTVWY";
int FRAME_LEN=DEFAULT_FRAME_LEN;	// 50 step 1
int SLIDE_STEP=10;//default value
int PHASED_START=0;//for Excel ouput - position of phased sequences relative to transcription start
int stdin_flag=0;
int show_help=0;
int detail_report=0;
int econom_report=0;
int special_report=0;
//int signal_report=0;
int user_max_level=2; //updated 2004
int user_min_level=0;

int count_method=1; //method - one sequence by default
int best_function=0;
int use_direct_search=0;
int use_inverse_search=0;
int use_symmetry_search=0;
int use_complementary_search=0;
int use_reverse=0;
int use_forward=0;
int whole_sequencies=0;
int NSIZE=4;
int ALPHABET_LEN=4;
char alphabet_pointer[30]="atgc";//=default_alphabet;
char alphabet_legend[30]="atgc";//default alphabet legend
static char main_name[64]="main_name.seq";
char second_name[64]="second_name.se2";
char res_name[64]="*.seq";
char rep_name[64]="*.rep";
char deffunc[MAX_ALPHABET_LEN+1]="tacg";
byte FR[MAX_ALPHABET_LEN]={1,0,3,2,4,5,6,7,8,9,10,11};//ATGC->TACG
byte FN[MAX_ALPHABET_LEN]={0,1,2,3,4,5,6,7,8,9,10,11};//ATGC->ATGC
char *abbrev[5]={"NW","D>","S<","C>","I<"};
char *hlp_text[]={"Complexity comparison - averaged profles",
		  "Use .exe filename.ext [switches]",
		  "Lempel-Ziv(LZ), Entropy(E), Entopy in Markov model(M), Wooton-Federhen(WF) and Linguistic complexity(LC)",
		  "These switches available",
		  "<Switches>",
		  "   -uX: set method: Z=Lempel-Ziv, E=Entropy, M=Markov model, W=Wooton-Federhen, L=Linguistic",
		  "   -aXX...: set alphabet (-aatgc default)",
//		  "   -b: find best function by other step",
		  "   -d -i -s -c: copying operations for Lempel-Ziv complexity (-d default)",
		  "      Letter d,i,s,c means these search methods:",
		  "      Direct, Inverse, Symmetry & Complementary search methods",
//		  "   -zXX...: set func (default -dtacg, no effect with -b)",
		  "   -fXX set frame size (50 default), use only with sliding window -mp",
		  "   -tXX set slide step (50 default), use only with -mp",
//		  "   -eXX set start position to output for phased sequences (relative to TS etc.)(0 default)",
		  "   -h or -?: show this text",
		  "   -i: get data from standart input",
//		  "//   -k(filename): machine configuration file",
		  "   -m(p|w): set count method",
		  "      These methods available:",
		  "      -mp count averaged complexity profile by one measure (Lempel-Ziv by default)",
		  "          for phased sequences",
		  "      -mw compare all 5 complexity methods (LZ,Entropy,EMarkov,LWF,LC)",
		  "   -rs  (special) show all averaged profile by all 5 methods for phased sequence set",
		  "   -y (alternative - complexity method number) -y0 = LZ, -y4 = Linguistic complexity)",
		  "   -x (Markov model order (order 1 - dinucleotides by default)",
		  "   -a_* : set pre-defined alphabet:",
		  "   -aatgc, -a_WS([at][gc]=WS), -a_RY ([ag][tc]=RY), -a_MK ([ac][gt]=MK) for DNA",
		  "   -a_p20 (ACDEFGHIKLMNPQRSTVWY), -a_p2 (Hydrophobicity) -a_p3charge, -a_p3surface for amino acids",
		  "Example: .exe filename.seq -d -f50 -t10",
		  "Output automaticallt will be written to 'filename'.out",
		  "Attention: input data must be in FASTA format",
		  ""};
struct key_profile{
		  char *key_name;	//имя ключа в командной строке
		  char *initial;	//явное задание в квадратных скобках
		  char *legend;		//подпись ветвей в графическом файле
		  char *remark;		//полная расшифровка терминов
		  };
//key_profile currprofile={NULL,NULL,NULL,NULL};

key_profile defined_keys[HOW_KEYS]={
		   {"atgc","atgc","ATGC",
		    "ATGC - standard DNA alphabet"},
		   {"_WS","[at][gc]","WS",
		   "Weak-Strong DNA alphabet: W=A/T, S=G/C"},
		   {"_RY","[ag][tc]","RY",
		   "Purine-Pyrimidine DNA alphabet: R=A/G, Y=T/C"},
		   {"_MK","[ac][gt]","MK",
		   "amino-keto DNA alphabet: M=A/C, Y=G/T"},
		   {"_p20","ACDEFGHIKLMNPQRSTVWY","ACDEFGHIKLMNPQRSTVWY",
		    "standard 20 amino acid alphabet"},
		   {"_p2","[ailmfpwv][rndcqeghksty]","-+",
		    "Hydrophobicity:(-)hydrophobic, (+)hydrophilic"},
		   {"_p3charge","[rhk][ancqgilmfpstwyv][de]","+0-",
		    "Amino acid charge:(+)base, (0)neutral, (-)acid"},
		   {"_p3surface","[rndqehk][acgpstwy][ilmfv]","+0-",
		    "Amino acid surface:(+)outer, (0)ambivalent, (-)inner"}
		   };
		//built-in profiles (in HTML)
//key_profile *profile=&currprofile;

//******************************
//int FASTAreport=1;
//int MaxNumberofSeqs=MaxNumberofSeqs;//was 1000
int MaxNumSteps=2000; //for splice sites //1100;
char FASTAtable_name[64];
char SignalDiscoveryFilename[64]="orlov__";
char SignalDiscoveryFilename2[64]="orlov2__";
int ComplexityMethod=0;
int LClongestword=MAX_FRAME_LEN;
char *MethodName[MAX_METHODS+1]={"Lempel-Ziv (CLZ)", "Entropy (CE)","Entropy high order (CM)","Wooton-Federhen (CWF)","Linguistic (CL)",""};
int MaxNumSeqOutputProfiles=20; //maximal number of output profiles
int debugregime=0;
int MarkovModelOrder=2;
/*
char *do_path(char *pre_path, char *name)
{
	int pre_length = strlen(pre_path);
	int length = strlen(pre_path) + strlen(name);
	char *path = (char*)malloc(length +2);

	strcpy(path, pre_path);
	if (pre_path[pre_length -1] != '/') {
		path[pre_length] = '/';
		pre_length++;
		length++;
	}
	strcpy(path +pre_length, name);
	path[length] = '\0';
	return path;
}
*/
FILE *std_msg=stdout;
s32bit *all_sequencies;
s32bit how_sequencies;
int seq_init() //allocate max size of sequence (file)
{int i;
 how_sequencies=0;
 all_sequencies=new s32bit[MAX_SEQUENCIES];
 if(all_sequencies!=NULL)
   {
   for(i=0;i<MAX_SEQUENCIES;i++)all_sequencies[i]=0;
   return 1;
   }else return 0;
}
s32bit add_sequence(s32bit item)//inner index for new sequence (start)
{
 if(how_sequencies<MAX_SEQUENCIES)
   {
   all_sequencies[how_sequencies++]=item;
   }
 return how_sequencies-1;
}
info determine_sequence(info item)//absolute address -> relative adress + number of sequence
{s32bit left,right,medium;
 left=0;
 right=how_sequencies-1;
 do
  {
  medium=(left+right)/2;
  if(item.from>=all_sequencies[medium])left=medium;
    else right=medium;
  }while((right-left)>1);
 if(all_sequencies[right]<=item.from)item.seq_number=right;
   else item.seq_number=left;
 item.from-=all_sequencies[item.seq_number];
 return item;
}
void seq_done()
{
 if(all_sequencies!=NULL)delete all_sequencies;
}
//******************************
//int handle_file(FILE *it);
int handle_key(char *key)
{int i;
// FILE *h;
  if(key[0]==SWITCH_CHAR)
   {
   switch(lowcase(key[1]))
	 {
//	 case 'a':for(i=0;(i<30)&&(key[i+2]);i++)alphabet_pointer[i]=key[i+2];
//		  alphabet_pointer[i]=0;
//		  break;
	 case 'a':i=0;
		  while((i<HOW_KEYS)&&(strcmp(key+2,defined_keys[i].key_name)))i++; //compare to standard alphabets ATGC, WS, etc.
		  if(i<HOW_KEYS)
		    {
		   // profile=&defined_keys[i];
		    strcpy(alphabet_pointer,defined_keys[i].initial);
		    strcpy(alphabet_legend,defined_keys[i].legend);
		    }
		    else
		    {
			for(i=0;(i<30)&&(key[i+2]);i++)alphabet_pointer[i]=key[i+2];
		    alphabet_pointer[i]=0;
		    }
		  break;

     case 'd':use_direct_search=1;
				      break;
     case 'i':use_inverse_search=1;
				      break;
	 case 's':use_symmetry_search=1;
				      break;
	 case 'c':use_complementary_search=1;
				      break;

	 case 'f':FRAME_LEN=atoi(key+2);
		  if(FRAME_LEN<MIN_FRAME_LEN) FRAME_LEN=MIN_FRAME_LEN;
		  if(FRAME_LEN>MAX_FRAME_LEN) FRAME_LEN=MAX_FRAME_LEN;
		  break;
	 case 't':SLIDE_STEP=atoi(key+2);
		  if(SLIDE_STEP<1) SLIDE_STEP=1;
		  if(SLIDE_STEP>MAX_FRAME_LEN) SLIDE_STEP=MAX_FRAME_LEN;
		  break;
//	 case 'e':PHASED_START=atoi(key+2);
//		  break;
     case 'z':if(strlen(key+2)>2) strcpy(deffunc,key+2); //at least 2 symbols for complementary
				else fprintf(std_msg,"Bad complementary function defined in command line. Use default: %s\n",deffunc);
		 break;
	 case 'h':
	 case '?':i=0;
		  show_help=1;
		  while(hlp_text[i][0])fprintf(std_msg,"%s\n",hlp_text[i++]);
		  exit(1);
		  break;
	 case 'm':switch(lowcase(key[2]))
			{
			case 'p':count_method=0;
				 break;
			case 'w':count_method=1;
				 break;
			default:  		 
				 fprintf(std_msg,"Bad data processing method defined in command line. Use default: %d\n",count_method);
				 break;
			}
		  break;
	 case 'u':switch(lowcase(key[2]))
			{
			case 'z':ComplexityMethod=0;//Lempel-Ziv
				break;
			case 'e':ComplexityMethod=1;//Entropy
				break;
			case 'm':ComplexityMethod=2;//Markov model (order>1)
				break;
			case 'w':ComplexityMethod=3;//Wooten-Federhen
				break;
			case 'l':ComplexityMethod=4;//Linguistic
				break;
			default:  		 
				fprintf(std_msg,"Bad complexity method defined in command line. Use default");
				break;
			}
		  break;
	 case 'r':
		 switch(lowcase(key[2]))       //detailed output may be with restrictions
			{
			case 'f':
					detail_report=1;
				 break;
			case 'e':
					econom_report=1;
				  break;
			case 's':
					special_report=1;
				 break;
			case 'd'://debug regime
					debugregime=1;
				  break;

			case 'n':
			default: 
				detail_report=0;
				econom_report=0;
				special_report=0;
				 break;
			} 
		 break;
//	 case 'o':strcpy(res_name,key+2);
//		  break;
	 case 'y': //define min (y) and max (x) for 'econom' output
			{ 
			user_min_level=atoi(key+2);
			if(user_min_level<0)user_min_level=0;
		//	if(user_min_level>MAX_FRAME_LEN)user_min_level=MAX_FRAME_LEN-1;
			ComplexityMethod=user_min_level;
			}
		  break;
	 case 'x': 
			{ 
			user_max_level=atoi(key+2);
		//	if(user_max_level>MAX_FRAME_LEN)user_max_level=MAX_FRAME_LEN;
			user_max_level++;
			if(user_max_level<2)user_max_level=2;//Markov order model
			MarkovModelOrder=user_max_level;
			}
		  break;
	 case 'k': //define min (y) and max (x) for 'econom' output
			{ 
			LClongestword=atoi(key+2);
			if(LClongestword<2)LClongestword=2;
			if(LClongestword>MAX_FRAME_LEN)LClongestword=MAX_FRAME_LEN;
		//	if(user_min_level>MAX_FRAME_LEN)user_min_level=MAX_FRAME_LEN-1;
			}
		  break;
	 case 'n':NSIZE=atoi(key+2);
		  if(NSIZE<1) NSIZE=1;
		  if(NSIZE>9) NSIZE=9;
		  break;
	 case 'w': 
			{ 
			MaxNumSeqOutputProfiles=atoi(key+2);
		//	if(user_max_level>MAX_FRAME_LEN)user_max_level=MAX_FRAME_LEN;
		//	user_max_level++;
			if(MaxNumSeqOutputProfiles<0)MaxNumSeqOutputProfiles=0;
			if(MaxNumSeqOutputProfiles>MAX_SEQUENCIES)MaxNumSeqOutputProfiles=MAX_SEQUENCIES;
			}
		  break;

	  default:return 2;
		  break;
	 }
   return 0;
   }else return 1;
}

//******************************
genfileinput *base;
int fill_sequence(byte *buff,s32bit *how)//read next (one) sequence to buffer from file
{int result=base->get_code();
 s32bit i=0;
 while((result!=END_OF_DATA)&&(result!=NEW_SEQUENCE)&&(i<SEQUENCE_LEN))
      {
      buff[i]=result;
      result=base->get_code();
      i++;
      }
 *how=i;
 return result;
}
char show_item_str[MAX_ALPHABET_LEN+2];
void show_item(info it,byte *buffer,int index)	//show one component of decomposition 
{s32bit i;
 fprintf(std_msg,"%i ([%i:%i],%i",index+1,it.seq_number+1,it.from+1,it.len);
 if(best_function)								//show for best variant decomposition 
   {
   for(i=0;i<ALPHABET_LEN;i++)show_item_str[i]='?';
   show_item_str[ALPHABET_LEN+1]=0;
   for(i=0;i<ALPHABET_LEN;i++)if(FR[i]!=255)show_item_str[FR[i]]=base->legend[i];
   if(it.method==3)show_item_str[ALPHABET_LEN]='>';
     else show_item_str[ALPHABET_LEN]='<';
   fprintf(std_msg,",%s,'",show_item_str);
   }else fprintf(std_msg,",%s,'",abbrev[it.method]);
 if(it.len)
   {
   for(i=0;i<it.len;i++)fprintf(std_msg,"%c",base->legend[buffer[i]]);
   }else fprintf(std_msg,"%c",base->legend[buffer[0]]);
 fprintf(std_msg,"')\n");
 //fflush(std_msg);
}
abstract_divizion *work;
//char tempfilename[64]="data.tmp";
info find_maximal_match(byte *buff,int size,int limitf,int limitr) //Heart !!!
{info h0={0,0,0},h1;
 int i;
 if(best_function)
   {
   if(use_forward)
     {
     h1=work->find_forward_match(buff,size,limitf,FR);
     if(h1.len>h0.len)h0=h1;
     }
   if(use_reverse)
     {
     h1=work->find_reverse_match(buff,size,limitr,FN);
     if(h1.len>h0.len)
       {
       for(i=0;i<MAX_ALPHABET_LEN;i++)FR[i]=FN[i];
       h0=h1;
       }
     }
   }else
   {
   if(use_direct_search)
     {
     h1=work->find_forward_match(buff,size,limitf,FN);
     if(h1.len>h0.len)h0=h1;
	   else if((h1.len==h0.len)&&(h1.from>h0.from))h0=h1; //select component more far the start
     }
   if(use_symmetry_search)
     {
     h1=work->find_reverse_match(buff,size,limitr,FN);
     if(h1.len>h0.len)h0=h1;
	    else if((h1.len==h0.len)&&(h1.from>h0.from))h0=h1;
     }
   if(use_complementary_search)
     {
     h1=work->find_forward_match(buff,size,limitf,FR);
     if(h1.len>h0.len)h0=h1;
	   else if((h1.len==h0.len)&&(h1.from>h0.from))h0=h1;
     }
   if(use_inverse_search)
     {
     h1=work->find_reverse_match(buff,size,limitr,FR);
     if(h1.len>h0.len)h0=h1;
	   else if((h1.len==h0.len)&&(h1.from>h0.from))h0=h1;
     }
   }
 if(!h0.len)h0.len++;
 return h0;
}
byte *buff0,*buff1;
void show_item_list(list_info* head,byte* universe,int start) //-- for profile for window
{int complexity=0;
 while(head!=NULL)
	{
	show_item(head->dt,universe+start,start);
	complexity++;
	start+=head->dt.len;
	head=head->next;
	}
 fprintf(std_msg,"Complexity is %i\n",complexity);
}
int ft_vector2=0;
list_info* new_decomposition(list_info *head,int* complexity,int number)//Nemytikova's algorithm
{int i=0,old_pnt=0,comple=0;
 list_info res;
 list_info *result=&res,*hlp,*old_work=head->next;
 ft_vector2++;
 byte *dt=buff0+ft_vector2;
 while(i<FRAME_LEN)
	{
	result->next=new list_info;
	result=result->next;
	if((i==old_pnt)&&(old_work->next!=NULL)&&(old_work->dt.from>=ft_vector2))
		  {
		  result->dt=old_work->dt;
		  }else
	result->dt=find_maximal_match((byte*)dt+i,FRAME_LEN-i,i+ft_vector2,i+ft_vector2);
	result->dt.seq_number=number;
	comple++;
	if(!result->dt.method)result->dt.from=i+ft_vector2;
	result->next=NULL;
	i+=result->dt.len;
	if(old_work->next!=NULL)
	  {
	  if(old_work->dt.len+old_pnt<=i)
		{
		old_pnt+=old_work->dt.len;
		old_work=old_work->next;
		}
	  }// end if(old_work->next!=NULL)
	}
 while(head!=NULL)
	{
	hlp=head->next;
	delete head;
	head=hlp;
	}
 *complexity=comple;
 return res.next;
}
int how_items=0;
float average=0;
float dispersion=0;
float summ=0;
float qsumm=0;
int max_item=0;
int min_item=MAX_FRAME_LEN;
float sigma_3=0;
void clean_statistics()
{
 how_items=0;
 average=0;
 dispersion=0;
 summ=0;
 qsumm=0;
 max_item=0;
 min_item=MAX_FRAME_LEN;
 sigma_3=0;
}
void collect_statistics(int it)
{
 if(it>max_item)max_item=it;
 if(it<min_item)min_item=it;
 summ+=(float)it;
 qsumm+=(float)(it*it);
 how_items++;
}
void show_statistics(FILE *rep) // statistics for sequence
{
 fprintf(rep,"Counted items:%i\n",how_items);
 if(how_items)
   {
   fprintf(rep,"Range:[%i;%i]\n",min_item,max_item);
   average=summ/(float)how_items;
   fprintf(rep,"Average:%f\n",average);
   if(how_items>1)
     {
	 dispersion=(qsumm-(summ*summ)/(float)how_items)/(float)(how_items-1);
	 fprintf(rep,"Dispersion:%f\n",dispersion);
	 sigma_3=dispersion*(float)2.5;
	 fprintf(rep,"3-sigma range:[%f;%f]\n",average-sigma_3,average+sigma_3);
     }   
   }
}

int calc_dispers_etc2(double *massiv, int n, double& resmin, double& resmax, double& resaverage, double& resdispers) //calculate dispersion, min, max, average
{
 double sum_mass=0;
 double qsum_mass=0;
 double aver_mass=0;
 double disp_mass=0;
 double tempi;
 double max_mass;
 double min_mass;
 int i;
 if(n<=0 || n>MaxNumberofSeqs) {printf("dispersion calculation error, n=%i\n",n);exit(1);}

 max_mass=massiv[0];
 min_mass=max_mass;
 for(i=0;i<n;i++)
	{
	tempi=massiv[i];
    if(tempi>max_mass)max_mass=tempi;
    if(tempi<min_mass)min_mass=tempi;
    sum_mass+=(double)(tempi);
    qsum_mass+=(double)(tempi*tempi);
	}
 aver_mass=sum_mass/(double)n;
 if(n>1)disp_mass=sqrt((qsum_mass-((sum_mass*sum_mass)/(double)n ) )/(double)(n-1) );
 else disp_mass=0;
	 resmax=max_mass;
resmin=min_mass;
resaverage=aver_mass;
resdispers=disp_mass;

return 0;
}


int compare_spec(const void *arg1, const void *arg2 );
int compare_spec(const void *arg1, const void *arg2 )
   {
   return ((int *)(arg1)-(int *)(arg2));
   }

//--add April
char ito_atgc(int i)
{char ch;
switch(i)
{
case 0: ch='a'; break;
case 1: ch='t'; break;
case 2: ch='g'; break;
case 3: ch='c'; break;
default:
	ch='X';
		break;
}
return(ch);
}

double estimate_entropy(byte *buff, int len, int ALPH_LEN)
{int count[256],i,countN=0;
double ent=0,tmp=0;

if(ALPH_LEN>256 || ALPH_LEN<1) {return(0);printf("(ent) bad ALPH_LEN=%i?\n",ALPH_LEN);}

if(len<1) {return(0);printf("(ent) bad length len=%i?\n",len);}
for(i=0;i<ALPH_LEN;i++) count[i]=0;//counts for nucleotides
for(i=0;i<len;i++)
{
 if((int)buff[i]<ALPH_LEN)count[(int)buff[i]]++;
 else countN++;//{printf("non ATGC?\n"); }
}
len-=countN;
 if(len<2) return(1);//many N, maximal entropy

for(i=0;i<ALPH_LEN;i++) 
{
tmp=(double)count[i];//+(double)(countN)/ALPH_LEN;
if( tmp>0)
  ent-=( tmp/(double)(len) )*log( tmp/(double)(len) );
}

//1./(log((double)ALPH_LEN)*ln1_2) )
//return (ent/(log( (double)(ALPH_LEN) )*ln1_2) );
//return (ent/(ln1_2) );
return (ent/(log( (double)(ALPH_LEN) ) ) );
}

double estimate_entropy2(byte *buff, int len, int ALPH_LEN)
{int count[400],i,len2;//Alphabet<16!
int countNN=0;
double entropy2=0,tmp;

if(ALPH_LEN>20 || ALPH_LEN<1) {return(0);printf("(ent2) bad ALPH_LEN=%i?\n",ALPH_LEN);}
if(len<1) {return(0);printf("(ent) bad length len=%i?\n",len);}
for(i=0;i<ALPH_LEN*ALPH_LEN;i++) count[i]=0;//counts for dinucleotides

for(i=0;i<len-1;i++)
{
 if((int)(buff[i]<ALPH_LEN)&&((int)buff[i+1]<ALPH_LEN))
	{count[((int)buff[i])*ALPH_LEN+((int)buff[i+1])]++;}
 else  countNN++;//{{printf("non ATGC?\n");}
}
len2=len-1-countNN;
 if(len2<2) return(1);//many NN, maximal entropy

for(i=0;i<ALPH_LEN*ALPH_LEN;i++) 
{
tmp=(double)count[i];
if( tmp>0)
 entropy2-=(tmp/(double)(len2))*log( tmp/(double)(len2));
}

//return (entropy2/(log( (double)(ALPH_LEN*ALPH_LEN)) *ln1_2 ) );
//return (entropy2/( ln1_2 ) );
return (entropy2/(log( (double)(ALPH_LEN*ALPH_LEN)) ) );
}

double estimate_entropy2more(byte *buff, int len, int ALPH_LEN, int order)
{int count[8000],i,j,len2;//Alphabet<20!
int countNN=0;
double entropy3=0,tmp;
int A_P_ORDER=1,A_P_ORDER_1=1,itmp=1;

if(ALPH_LEN>20 || ALPH_LEN<1) {return(0);printf("(ent2) bad ALPH_LEN=%i?\n",ALPH_LEN);}
if(len<=order) {return(0);printf("(ent 3) bad length %i < order len=%i?\n",len,order);}
if(order<2){return(0);printf("(ent 3) bad order=%i?\n",order);}

for(i=0;i<order-1;i++) A_P_ORDER_1*=ALPH_LEN;

A_P_ORDER=A_P_ORDER_1*ALPH_LEN;

for(i=0;i<A_P_ORDER;i++) count[i]=0;//counts for dinucleotides


for(i=0;i<=len-order;i++)
{
//index=(index*ALPHABET_LEN+currcode)%how_context;
 if((int)buff[i+order]<ALPH_LEN)
 {
 itmp=(int)buff[i];
  for(j=1;j<order;j++) //first oligonucleotide
    itmp=itmp*ALPH_LEN+((int)buff[i+j]);

 count[itmp]++;//((int)buff[i])*ALPH_LEN+((int)buff[i+1])]++;}
 }
 else  countNN++;//{printf("non ATGC? = %i\n",(int)(buff[i+order]));}
}
len2=len+1-order-countNN;

if(len2<2) return(1);//many NN, maximal entropy

for(i=0;i<A_P_ORDER;i++) 
{
tmp=(double)count[i];
if( tmp>0)
 entropy3-=(tmp/(double)(len2))*log( tmp/(double)(len2) );
}

return (entropy3/log((double)(A_P_ORDER)) );
}


int calc_freq(byte *buff, int len, int ALPH_LEN, int *count)
{int i;//,len2;//Alphabet<16!count[256],
int countN=0;
double entropy2=0;//,tmp;

if(ALPH_LEN>256 || ALPH_LEN<1) {return(0);printf("(calc.freq) bad ALPH_LEN=%i?\n",ALPH_LEN);}
if(len<1) {return(0);printf("(calc.freq) bad length len=%i?\n",len);}

for(i=0;i<ALPH_LEN;i++) count[i]=0;//counts for nucleotides

for(i=0;i<len;i++)
{
 if((int)buff[i]<ALPH_LEN)count[(int)buff[i]]++;
 else countN++;//{printf("non ATGC?\n"); }
}

return (countN);
}


int calculateFASTAprofile(FILE *f0,FILE *rep) //main profile calculation
{int result,how,from,size,number=0,i=0,j=0,NN;
 list_info hdr={{0,0,0,0},NULL};
 list_info* resu=&hdr;
 info h0;

int FASTAcurrstep=0;//position in profile
int FASTAmaxstep=1; //maximal profile length(if different sequences' lengths)
double **FASTAtable;
double **MethodTable;//for complexity methods comparison
double *FASTAmeanpos;
double *FASTAminpos,*FASTAmaxpos; int *FASTAposnon0;
int numallseq=0;
int MaxNumSeqsForTable=1000;//make shorter table for output to Excel
int FASTAmintablelength=MaxNumSteps;
double complexity;//, complexityLn, complexityLZ;
double words9=0, Entropy=0, Entropy2=0, Entropy9=0,cmet=0,cWooten=0.,cLC=0,cLCbigK=0;
int kk; int *countfreq; //int A_POWER_NSIZE;
int *counteq,*countmax;
//FILE *fastaout;
//FILE *fastaout2;
double LCmax=0,LCmaxbigK=0;
countfreq=new int[ALPHABET_LEN];
countmax=new int[FRAME_LEN+1];
counteq=new int[FRAME_LEN+1];

//countmax=new int[FRAME_LEN+1];

countmax[0]=0;
i=0;how=ALPHABET_LEN;
while(how<FRAME_LEN-i+1) 
	{
	countmax[i+1]=how;
    how*=ALPHABET_LEN;
	i++;
    }
for(kk=i+1;kk<FRAME_LEN;kk++)
	countmax[kk]=FRAME_LEN-kk+1;

//for(i=1;i<=FRAME_LEN;i++)
//	LCmaxbig+=(double)(countmax[i]);
if(LClongestword>FRAME_LEN)LClongestword=FRAME_LEN;
if(LClongestword<NSIZE)LClongestword=NSIZE;
for(i=1;i<=LClongestword;i++)
	LCmaxbigK+=(double)(countmax[i]);

if(MarkovModelOrder>NSIZE)MarkovModelOrder=NSIZE;


//ComplexityMethod=user_min_level;

if(MaxNumSeqsForTable>MAX_SEQUENCIES)MaxNumSeqsForTable=MAX_SEQUENCIES;

FASTAtable=new double*[MAX_SEQUENCIES+1];
		FASTAmeanpos=new double[MaxNumSteps+1];
		FASTAminpos=new double[MaxNumSteps+1];
		FASTAmaxpos=new double[MaxNumSteps+1];
		FASTAposnon0=new int[MaxNumSteps+1];
		countfreq=new int[ALPHABET_LEN];
		for (i=0;i<MaxNumSteps;i++)
			{
			FASTAmeanpos[i]=0.;//initialize, maybe not need due to (new) operator
			FASTAminpos[i]=-1;
			FASTAmaxpos[i]=0;
			FASTAposnon0[i]=0;
			}
if(special_report)
		{
		MethodTable=new double*[MAX_METHODS];//[0][FASTAcurrstep]
		if(MethodTable!=NULL)
			{
			for (i=0;i<MAX_METHODS;i++)
				MethodTable[i]=new double[MaxNumSteps+1];
			}
		}
if(FASTAtable!=NULL)
	{ 
//	for (i=0;i<MaxNumberofSeqs+1;i++)
	for (i=0;i<MaxNumSeqsForTable;i++)
		{
		FASTAtable[i]=new double[MaxNumSteps+1];
		if(FASTAtable[i]==NULL){fprintf(std_msg,"Not enough memory for table (%i x %i) for profiles output\n",MaxNumberofSeqs,MaxNumSteps);exit(1);}
		for (j=0;j<MaxNumSteps;j++) FASTAtable[i][j]=0;
		}
	}
else {fprintf(std_msg,"Not enough memory for tables for simultaneous profiles output\n");exit(1);}
//}

i=0;
j=0;
how=0;
	   if(SLIDE_STEP<MIN_SLIDE_STEP)SLIDE_STEP=MIN_SLIDE_STEP;
//fprintf(std_msg,"Report(signal=%i)Frame length:%i\nSlide step:%i\n",signal_report,FRAME_LEN,SLIDE_STEP);

 base->set_file(f0);
//--------------- main cycle
 do{
   result=fill_sequence(buff0,&how); //find and set next sequence in the file (in FASTA format)
   if(!how)continue;
   if((how>0) && (how<FRAME_LEN)) continue;
//if(debugregime)
//	{
//	fprintf(std_msg,"Debug!!! how=%i FRAME_LEN=%i SLIDE_STEP=%i\n",how,FRAME_LEN,SLIDE_STEP);
//	}
 
 if(FASTAmintablelength>(int)((double)(how-FRAME_LEN)/(double)SLIDE_STEP)) 
 {
	FASTAmintablelength=(int)((double)(how-FRAME_LEN)/(double)SLIDE_STEP);
	if(number==0) fprintf(std_msg,"Standard Output to table %i lines (real seq.length=%i)\n",FASTAmintablelength,FASTAmintablelength*SLIDE_STEP+FRAME_LEN);
    if(number>0)  fprintf(std_msg,"Corrected table size %i lines (real seq.length=%i)\n",FASTAmintablelength,FASTAmintablelength*SLIDE_STEP+FRAME_LEN);
 }

 if(FASTAmaxstep<FASTAcurrstep) FASTAmaxstep=FASTAcurrstep;
	FASTAcurrstep=0;

   from=0;
     do{
       if((from+FRAME_LEN)<=how)size=FRAME_LEN;
	 else size=how-from;

       work->clear();
       work->set_universe(buff0+from,size);
       int limit=0;
       complexity=0;
	   words9=0; Entropy=1; Entropy2=1; Entropy9=1;
	   i=0;

       do{
	 if(i>=limit)
	   {
	   h0=find_maximal_match(buff0+from+i,size-i,limit,limit);
	   h0.seq_number=number;
	   limit+=h0.len;
	   complexity++;
	   }
	 work->add_vector(buff0+from+i,i,size-i);
	 i++;
	 }while(i<size);

     cLCbigK=(double)work->linguistic_bigK(buff0+from, size, counteq, LClongestword);
//     cLCbig=(double)work->linguistic_big(buff0+from, size, counteq);
     for(kk=NSIZE+1;kk<=LClongestword;kk++)
		{
		 if((countmax[kk]-counteq[kk])<0) {counteq[kk]=countmax[kk];printf("bug %i>%i\n",counteq[kk],countmax[kk]);}//flagdebug=1;
			if(countmax[kk]==1) counteq[kk]=0;//last number 1 - n in n
	    cLCbigK+=(countmax[kk]-counteq[kk]);
		}
	 cLCbigK/=LCmaxbigK;

//-------- Entropy
	 Entropy=estimate_entropy(buff0+from, size, ALPHABET_LEN);
	 Entropy9=estimate_entropy2more(buff0+from, size, ALPHABET_LEN, MarkovModelOrder);
//-------- Wooten 
	 NN=0;
NN=calc_freq(buff0+from, size, ALPHABET_LEN, countfreq);
if(size>NN) //not all symbols N
{
cWooten=logfak(size-NN);
for(kk=0;kk<ALPHABET_LEN;kk++)
 cWooten-=logfak(countfreq[kk]);
	 //K=1/L logN (L! /  ni!)
 cWooten/=(double)(size-NN)*( log((double)ALPHABET_LEN) );
}
else
 cWooten=1;

switch(ComplexityMethod)
{
	case 0: 
		//cmet=complexity/FRAME_LEN; //complexity=complexity - number of operations
        //cmet=(complexity)/(log((double)ALPHABET_LEN)*log((double)FRAME_LEN ));//1-((double)FRAME_LEN / complexity)/(log((double)ALPHABET_LEN)*log((double)FRAME_LEN ));
		cmet=complexity/((double)(FRAME_LEN-NN) );
		break;
	case 1: 
		cmet=Entropy;
		break;
	case 2: 
		cmet=Entropy9;
		break;
	case 3: 
		cmet=cWooten;
		break;
	case 4: 
		cmet=cLCbigK;
		break;
	default:
		cmet=complexity;
		break;
}

       from+=SLIDE_STEP;

	   if(special_report) 
	   {
		MethodTable[0][FASTAcurrstep]+= complexity/(double)(FRAME_LEN-NN) ;
		MethodTable[1][FASTAcurrstep]+=Entropy;
		MethodTable[2][FASTAcurrstep]+=Entropy9;
		MethodTable[3][FASTAcurrstep]+=cWooten;
		MethodTable[4][FASTAcurrstep]+=cLCbigK;
	   }

	if(number<MaxNumSeqsForTable) FASTAtable[number][FASTAcurrstep]=cmet;
	if(cmet>0)
		{
	    FASTAmeanpos[FASTAcurrstep]+=(double)cmet; //first window 
		FASTAposnon0[FASTAcurrstep]++;
	    if(FASTAminpos[FASTAcurrstep]>cmet)FASTAminpos[FASTAcurrstep]=cmet;
	    if(FASTAmaxpos[FASTAcurrstep]<cmet)FASTAmaxpos[FASTAcurrstep]=cmet;
	    if(FASTAposnon0[FASTAcurrstep]==1) //inializing profile parameters from first (NON-ZERO) sequence
			{
			FASTAminpos[FASTAcurrstep]=cmet;
			FASTAmaxpos[FASTAcurrstep]=cmet;
			}
		}
    if(cmet==0) 
		{
		printf("cmet=0,seq. num=%i position*SLIDE_STEP=%i, FRAME_LEN=%i\n",number,FASTAcurrstep*SLIDE_STEP,FRAME_LEN);//old check - debug
		if(FASTAmintablelength>(int)((FASTAcurrstep-FRAME_LEN)/SLIDE_STEP)) 
			{
			FASTAmintablelength=(int)((FASTAcurrstep-FRAME_LEN)/SLIDE_STEP);//mark minimal seq.length in the sample
			fprintf(std_msg,"Corrected for shorter seq-s(non-ATGC probably) table %i lines (real seq.length=%i)\n",FASTAmintablelength,FASTAmintablelength*SLIDE_STEP+FRAME_LEN);
			}
		}

    if(FASTAcurrstep<MaxNumSteps)	FASTAcurrstep++;

       }while((from+size)<how);//here was < 

   if(FASTAmaxstep<FASTAcurrstep) FASTAmaxstep=FASTAcurrstep;
   number++;
//   if(number==MaxNumSeqsForTable) {printf("only %i sequencies for output\n",number);}//exit(1);}
   numallseq++; //non-limited number of sequences - to calculate mean values only
   }while(result!=END_OF_DATA && number<MaxNumberofSeqs);

//if(FASTAreport && (!signal_report))//march 2003 -updated July 2003-updated - Aug+signals to Discovery
//{ 
FILE *fastaout;
//FILE *fastaout3;
//double tempstat;
double tempmin,tempmax;
double tempaverage, tempdispers=0;
double *tempmassiv;//, *histomassiv;
int numbernon0=0;
tempmassiv=new double[MAX_SEQUENCIES];
//histomassiv=new double[50];//histogram of complexity distribution in point for Excel

//if( (fastaout=fopen(do_path(task_path,FASTAtable_name),"w"))==NULL) {fprintf(std_msg,"file FASTA profiles output is not opened\n");exit(1);}
if( (fastaout=fopen(FASTAtable_name,"w"))==NULL) {fprintf(std_msg,"file FASTA profiles output is not opened\n");exit(1);}
//if( (fastaout3=fopen(do_path(task_path,"histo2.out"),"w"))==NULL) {fprintf(std_msg,"file FASTA profiles2 output is not opened\n");exit(1);}

if(number>MaxNumSeqsForTable) number=MaxNumSeqsForTable; //correction for table output
printf("OK Seq.number=%i,(MaxNumSeqsForTable=%i),numallseq=%i,MaxNumberofSeqs=%i\n",number,MaxNumSeqsForTable,numallseq,MaxNumberofSeqs);

if(special_report)
	{
//	FILE *fastaout2;
//	if( (fastaout2=fopen(do_path(task_path,"comparemeanprofiles5.out"),"w"))==NULL) {fprintf(std_msg,"file FASTA profiles2 output is not opened\n");exit(1);}
	if(number<1) fprintf(fastaout,"sorry, no sequences\n");
	else
		{
		fprintf(fastaout,"Averaged profiles for set of sequences by all complexity methods\n");
		fprintf(fastaout,"Pos.\tCLZ\tCE\tCM(%i)\tCWF\tCL(%i)\n",MarkovModelOrder,LClongestword);
	//	fprintf(fastaout,"Pos.\tLZ\tEnt.\tEnt.h\tLWF\tLC\n");//again average 
		for(i=0;i<FASTAmaxstep;i++)
			{
			fprintf(fastaout,"%i",(int)((FRAME_LEN/2)+i*SLIDE_STEP+PHASED_START));
			for(j=0;j<MAX_METHODS;j++)
				{
				if(FASTAposnon0[i]>0) MethodTable[j][i]/=(double)(FASTAposnon0[i]);
				fprintf(fastaout,"\t%.3f",MethodTable[j][i]);
				}
			fprintf(fastaout,"\n");
			}
		}
	}


if(!special_report)
{
fprintf(fastaout,"Total number of sequences analysed by method %s is %i\n",MethodName[ComplexityMethod],number);//again average 
fprintf(fastaout,"Averaged complexity profile for (only first %i seq. in columns)\n",(number<MaxNumSeqOutputProfiles)?number:MaxNumSeqOutputProfiles);//again average 
fprintf(fastaout,"Pos.\tMean\tMin\tMax\tSigma");//again average 
for(i=0;i<((number<MaxNumSeqOutputProfiles)?number:MaxNumSeqOutputProfiles);i++)
 fprintf(fastaout,"\tSeq%i",i+1);//seq.numbers in table columns

fprintf(fastaout,"\n");//seq.numbers in table columns

for(i=0;i<FASTAmaxstep;i++)
	{
//	fprintf(fastaout,"%i",(FRAME_LEN/2)+i*SLIDE_STEP);
	if(FASTAposnon0[i]<1) {FASTAposnon0[i]=1;printf("bug? FASTAposnon0[i]<1,FASTAcurrstep=%i\n",i);}//exit(1);}
	if(FASTAposnon0[i]>numallseq) {printf("bug? FASTAposnon0[i]>numallseq,i=%i\n",i);exit(1);}
	if(FASTAposnon0[i]!=0)	 FASTAmeanpos[i]/=(double)(FASTAposnon0[i]); 
	fprintf(fastaout,"%i\t%.3f\t%.3f\t%.3f",(int)((FRAME_LEN/2)+i*SLIDE_STEP+PHASED_START),FASTAmeanpos[i],FASTAminpos[i],FASTAmaxpos[i]);   

	numbernon0=0;
	for(j=0;j<number;j++)
		{
		if(FASTAtable[j][i]>0) 
			{
			tempmassiv[numbernon0]=FASTAtable[j][i];
			numbernon0++;
			}
		}
if(numbernon0>0)
	calc_dispers_etc2(tempmassiv, numbernon0, tempmin, tempmax, tempaverage, tempdispers); //calculate dispersion, min, max, average
fprintf(fastaout,"\t%.2f",tempdispers);   //dispersion in the position

if(MaxNumSeqOutputProfiles>number)MaxNumSeqOutputProfiles=number;
	
 for(j=0;j<MaxNumSeqOutputProfiles;j++)//<number - here is less for Excel output
		fprintf(fastaout,"\t%.3f",FASTAtable[j][i]);

	fprintf(fastaout,"\n");//average profile
	}
//else	fprintf(fastaout,"\tshort sequences or no data\n");//average profile


}
fclose(fastaout);
//}  

 return 0;
}


int calculateLCW(FILE *f0,FILE *rep) //main profile calculation
{int result,how,from,size,number=0,i=0,j=0,NN=0,iLZ=0;
 list_info hdr={{0,0,0,0},NULL};
 list_info* resu=&hdr;
 info h0;
int FASTAcurrstep=0;//position in profile
int FASTAmaxstep=1; //maximal profile length(if different sequences' lengths)
double CompComp[MAX_METHODS];
int numallseq=0;
int MaxNumSeqsForTable=1000;//make shorter table for output to Excel
//int FASTAmintablelength=MaxNumSteps;
double complexity; //, complexityLn, complexityLZ;
double words9=0, Entropy=0, Entropy2=0, Entropy9=0,cmet=0,cWooten=0,cLC=0,cLCbig=0,cLCbigK=0;
int kk,k2; 
int countfreq[MAX_ALPHABET_LEN];// int A_POWER_NSIZE;
int *counteq,*countmax;
int flagdebug=0;
FILE *fastaout;
//FILE *fastaout2;
double LCmaxbig=0,LCmaxbigK=0;//LCmax=0,

int *dotlist, *dotfromlist;//temp
dotlist=new int[FRAME_LEN+1];
dotfromlist=new int[FRAME_LEN+1];
for(kk=i;kk<FRAME_LEN;kk++){dotlist[kk]=0;dotfromlist[kk]=0;}

countmax=new int[FRAME_LEN+1];
counteq=new int[FRAME_LEN+1];


//for(kk=i+1;kk<FRAME_LEN;kk++)
//	countmax[kk]=FRAME_LEN-kk+1;

//for(i=1;i<=FRAME_LEN;i++)
//	LCmaxbig+=(double)(countmax[i]);
//for(i=1;i<=LClongestword;i++)
//	LCmaxbigK+=(double)(countmax[i]);

//ComplexityMethod=user_min_level;
//if( (fastaout2=fopen(do_path(task_path,"comparemethods.out"),"w"))==NULL) {fprintf(std_msg,"file FASTA profiles2 output is not opened\n");exit(1);}
//if( (fastaout=fopen(do_path(task_path,FASTAtable_name),"w"))==NULL) {fprintf(std_msg,"file FASTA profiles output is not opened\n");exit(1);}
if( (fastaout=fopen(FASTAtable_name,"w"))==NULL) {fprintf(std_msg,"file FASTA profiles output is not opened\n");exit(1);}

if(MarkovModelOrder>NSIZE)MarkovModelOrder=NSIZE;

if(debugregime)
fprintf(fastaout,"Complexity estimations. Method comparison.\nALPHABET LENGTH=%i, l-gram size=%i\n",ALPHABET_LEN,NSIZE);

fprintf(fastaout,"Alphabet is ");
for(i=0;i<ALPHABET_LEN;i++) fprintf(fastaout,"%c",base->legend[i]);
if(!strcmp(alphabet_legend,defined_keys[0].legend)) //show 'N'; only for ATGC alphabet
	fprintf(fastaout," and undefined symbol=(%c)\n",base->legend[ALPHABET_LEN]);
else
	fprintf(fastaout,"\n");


//i=0;
j=0;

  if(SLIDE_STEP<MIN_SLIDE_STEP)SLIDE_STEP=MIN_SLIDE_STEP;

if(debugregime)
  fprintf(std_msg,"Frame length:%i\nSlide step:%i\n",FRAME_LEN,SLIDE_STEP);

 base->set_file(f0);
//--------------- main cycle
 do{
   result=fill_sequence(buff0,&how); //find and set next sequence in the file (in FASTA format)
   if(!how)continue;

if(how<FRAME_LEN)FRAME_LEN=how;//case of short sequence
if(LClongestword>how)LClongestword=how;
if(LClongestword>FRAME_LEN)LClongestword=FRAME_LEN;
if(LClongestword<1)LClongestword=1;
//if(NSIZE>FRAME_LEN)NSIZE=FRAME_LEN-1;
//if(NSIZE<1)NSIZE=1;

countmax[0]=0;
i=0;j=ALPHABET_LEN;
while(j<FRAME_LEN-i+1) 
	{
	countmax[i+1]=j;
    j*=ALPHABET_LEN;
	i++;
    }
for(kk=i+1;kk<=FRAME_LEN;kk++)
	countmax[kk]=FRAME_LEN-kk+1;
for(i=1;i<=LClongestword;i++)
	LCmaxbigK+=(double)(countmax[i]);


fprintf(fastaout,"Frame length:%i\tSlide step:%i\n",FRAME_LEN,SLIDE_STEP);
if(!econom_report)
	{
	fprintf(fastaout,"Complexity profiles. One sequence by length %i\n",how);
	fprintf(fastaout,"Pos.\tCLZ\tCE\tCM(%i)\tCWF\tCL(%i)\t%%NN\n",MarkovModelOrder,LClongestword);
	}
else
	{
	fprintf(fastaout,"Complexity profile. One sequence by length %i\n",how);
	fprintf(fastaout,"Pos.\t%s",MethodName[ComplexityMethod]);
	switch(ComplexityMethod)
		{
		case 2: 
		fprintf(fastaout," m=%i\n",MarkovModelOrder);
		break;
		case 4: 
		fprintf(fastaout," k=%i\n",LClongestword);
		break;
		default:
		fprintf(fastaout,"\n");
		break;
		}
//char *MethodName[MAX_METHODS+1]={"Lempel-Ziv (CLZ)", "Entropy (CE)","Entropy high order (CM)","Wooton-Federhen (CWF)","Linguistic (CL)",""};
	}


   from=0;

   do{
       if((from+FRAME_LEN)<=how)size=FRAME_LEN;
	 else size=how-from;
//Lempel-Ziv ---------------------------
       work->clear();
       work->set_universe(buff0+from,size);
       int limit=0;
       complexity=0;
//	   complexityLn=0;//	   complexityLZ=0;//	   words9=0; Entropy=1; Entropy2=1; Entropy9=1;
	   iLZ=0;

NN=0;
NN=calc_freq(buff0+from, size, ALPHABET_LEN, countfreq);
if(size>NN) //check for N symbols, if all N - skip
{

       do{
	 if(iLZ>=limit)
	   {
	   h0=find_maximal_match(buff0+from+iLZ,size-iLZ,limit,limit);
	   h0.seq_number=number;
	   limit+=h0.len;

	   //debugging - for test outpur
			dotlist[(int)complexity]=(int)(h0.len);
			dotfromlist[(int)complexity]=(int)(h0.from);
	   complexity++;
		}
	 work->add_vector(buff0+from+iLZ,iLZ,size-iLZ);
	 iLZ++;
	 }while(iLZ<size);

	if(size<=0) printf("bug size=%i\n",size);
	
	 //words9=work->count_low_level(FRAME_LEN,&Entropy9);
//     cLC=(double)work->linguistic(NSIZE);
//     cLC/=LCmax;
/*     cLCbig=(double)work->linguistic_big(buff0+from, size, counteq);
     for(kk=NSIZE+1;kk<=size;kk++)
		{
			if((countmax[kk]-counteq[kk])<0) counteq[j]=countmax[j];//flagdebug=1;
	    cLCbig+=(countmax[kk]-counteq[kk]);
		}
	 cLCbig/=LCmaxbig;
*/
     cLCbigK=(double)work->linguistic_bigK(buff0+from, size, counteq, LClongestword);

if((from==0)&&debugregime)
{
	printf("--- nsize=%i LClongestword=%i LCmaxbigK=%i----(FRAME_LEN=%i)-------------------\n",NSIZE,LClongestword,(int)LCmaxbigK,FRAME_LEN);
 k2=0;
	for(j=0;j<(int)complexity;j++)
	{
	for(kk=0;kk<dotlist[j];kk++) if((k2+kk)%10) printf("-"); else printf("|"); 
	printf("*");
	k2+=dotlist[j];
	}
	printf("\n");
	k2=0;
	for(j=0;j<(int)complexity;j++)
	{
	for(kk=0;kk<dotlist[j];kk++) printf("%c",base->legend[((int)buff0[from+k2+kk] )] );
	printf("*");
	k2+=dotlist[j];
	}
	printf("\n");
 k2=0;
	for(j=0;j<(int)complexity;j++)
	{
	for(kk=0;kk<(int)(dotlist[j]-1);kk++) {printf(" "); k2++;}
	printf("%2i",dotfromlist[j]);k2+=2;
	}
	printf("\n");
	 
	printf("FRAME_LEN=%i cLCbigK=%i LCMaxbigK(all)=%i\n",FRAME_LEN,(int)cLCbigK,(int)LCmaxbigK);
     for(j=0;j<=size;j++)
	   printf("%i\t",countmax[j]);
     printf("\n");
     for(j=0;j<=size;j++)
	   printf("%i\t",counteq[j]);
     printf("\n");
		kk=0;
	if(NSIZE<=LClongestword)
		{
     for(j=0;j<=NSIZE;j++)
		{printf("%i\t",counteq[j]);kk+=counteq[j];}
     for(j=NSIZE+1;j<=LClongestword;j++)
		{
		 if(countmax[j]==1) counteq[j]=0;// to substruct
		 printf("%i\t",countmax[j]-counteq[j]);kk+=countmax[j]-counteq[j];
		}
		}
		else
		{
		for(j=0;j<=LClongestword;j++)
			{printf("%i\t",counteq[j]);kk+=counteq[j];}
		}

     printf("\nkk=%i,LCmaxbigK=%i, kk/LCmaxbigK=%.3f\n",kk,(int)LCmaxbigK,(double)kk/(double)(LCmaxbigK));
     printf("\n");
	 flagdebug=0;
}
     for(kk=NSIZE+1;kk<=LClongestword;kk++)
		{
		 if((countmax[kk]-counteq[kk])<0) {counteq[kk]=countmax[kk];printf("bug %i>%i\n",counteq[kk],countmax[kk]);}//flagdebug=1;
			if(countmax[kk]==1) counteq[kk]=0;//last number 1 - n in n
	    cLCbigK+=(countmax[kk]-counteq[kk]);
		}
	 cLCbigK/=LCmaxbigK;


//-------- Entropy
	 Entropy=estimate_entropy(buff0+from, size, ALPHABET_LEN);
//     Entropy2=estimate_entropy2(buff0+from, size, ALPHABET_LEN);
	 Entropy9=estimate_entropy2more(buff0+from, size, ALPHABET_LEN, MarkovModelOrder);
//	 double estimate_entropy2more(byte *buff, int len, int ALPH_LEN, int order)

//-------- Wooten 
if(size>NN) //not all symbols N
{
//cWooten=(double)(logfak((int)(size-NN)) );
	cWooten=logfak(size-NN);
for(kk=0;kk<ALPHABET_LEN;kk++)
 cWooten-=logfak(countfreq[kk]);
	 //K=1/L logN (L! /  ni!)
 cWooten/=(double)(size-NN)*( log((double)ALPHABET_LEN) );
}
else
 cWooten=1;



//    case 0:    	 //complexity=complexity - number of operations
//CompComp[0]=(complexity)/(log((double)ALPHABET_LEN)*log((double)FRAME_LEN ));
//CompComp[0]=1-((double)FRAME_LEN / complexity)/(log((double)ALPHABET_LEN)*log((double)FRAME_LEN ));
 //CompComp[0]=((complexity/log((double)ALPHABET_LEN)*log((double)FRAME_LEN ) ))/((double)FRAME_LEN );

CompComp[0]=(complexity-(double)NN)/((double)(size-NN) );
CompComp[1]=Entropy;
CompComp[2]=Entropy9;
CompComp[3]=cWooten;
//CompComp[4]=cLCbig;
CompComp[4]=cLCbigK;
}
else //no symbols, only 'N', all complexity = 1
{
CompComp[0]=1; 
CompComp[1]=1;
CompComp[2]=1;
CompComp[3]=1;
CompComp[4]=1;
}
	fprintf(fastaout,"%i",(FRAME_LEN/2)+from);
if(!econom_report)
	{	
	for(kk=0;kk<MAX_METHODS;kk++)//<number - here is less for Excel output
		fprintf(fastaout,"\t%.3f",CompComp[kk]);
		//fprintf(fastaout,"\t%.3f",FASTAtable[kk][FASTAcurrstep]);
	fprintf(fastaout,"\t%.3f",(double)NN/(double)(FRAME_LEN));
	}
else
{
	fprintf(fastaout,"\t%.3f",CompComp[ComplexityMethod]);
	fprintf(fastaout,"\t%.3f",(double)NN/(double)(FRAME_LEN));
}
	fprintf(fastaout,"\n");

    from+=SLIDE_STEP;
//    if(FASTAcurrstep<MaxNumSteps)	FASTAcurrstep++;
       }
		while((from+size)<how);//here was < 

//   FASTAmaxstep=FASTAcurrstep;
   number++;
   }
   while(result!=END_OF_DATA && number<1);//only one sequence here

  
fclose(fastaout);
return 0;
}





//------------------------------------------------------------------------


/*
int calculate_vector(FILE *f0,FILE *rep) //for -mv (4)!*2=48 variants for each sequence in file
{int result,how,j,complexity,limit;
 unsigned long i,total=factorial(ALPHABET_LEN)*2;
 exact *functor=new exact(ALPHABET_LEN);
 byte *curr_func;
 info h0;
 base->set_file(f0);
 do{
   result=fill_sequence(buff0,&how);
   if(!how)continue;
   work->set_universe(buff0,how);
   functor->reset();
   fprintf(rep,"(");
   for(i=0;i<total;i++)
      {
      if(!(i&1))curr_func=functor->get_next();
	  ///erase it string !
	  for(int kkk=0;kkk<ALPHABET_LEN;kkk++)fprintf(rep,"%c",base->legend[curr_func[kkk]]);
	  fprintf(rep," ");
      work->clear();
      j=0;
      complexity=0;
      limit=0;
      do{
	if(j>=limit)
	  {
	  if(!(i&1))h0=work->find_forward_match(buff0+j,how-j,limit,curr_func);//forward - "chetnii", from 0
	    else h0=work->find_reverse_match(buff0+j,how-j,limit,curr_func);
	  if(detail_report)show_item(h0,buff0+j,j);
	  if(h0.len)limit+=h0.len;
	    else limit++;
	  //fprintf(rep,"%i\n",h0.len);
	  complexity++;
	  }
	work->add_vector(buff0+j,j,how-j);
	j++;
	}while(j<how);
      fprintf(rep,"%i ",complexity); //one of 48 variants
      }
   fprintf(rep,")\n");
   }while(result!=END_OF_DATA);
 delete functor;
 return 0;
}
*/
int set_complementary_function() //user-defined complementary function
{byte *new_fr=new byte[ALPHABET_LEN];
 int res=1,i,r;
 if(strlen(deffunc)==(unsigned)ALPHABET_LEN)
   {
   for(i=0;i<ALPHABET_LEN;i++)new_fr[i]=(byte)i; //fill 0,1,2...        length alphabet may change
   for(i=0;i<ALPHABET_LEN;i++)
      {
      r=base->get_value(deffunc[i]);//take symbol from defined string
      if(r>=0)new_fr[r]=i;//-1 not used, 
	else                  //exit from subroutine
	{
	res=0;
	break;
	}
      }
   }else res=0;
 if(res)for(i=0;i<ALPHABET_LEN;i++)FR[i]=new_fr[i];
 delete new_fr;
 return res;
}


int main(int argc,char *argv[])//int argco,char *argvo[]) //------- console variant !!!
{
	int i,j;
int flagname=0;
 FILE *fin=NULL,*fsec=NULL,*fout=stdout;

// int argc;
 //char** argv;
/*
 task_path = argvo[1];
 argc = argco -1;
 argv = (char**)malloc(argc*sizeof(char*));
 argv[0] = argvo[0];
 for (i=1; i<argc; i++) {
	         argv[i] = argvo[i+1];
 }
  */                                                                                                 
                                                                                                   
 

 if(argc>=2)
   {
	// analyse command line parameters
    for(i=1;i<argc;i++) //if(handle_key(argv[i])==1)strcpy(main_name,argv[i]);
	{
      switch(handle_key(argv[i]))	//
	    {
	    case 0:	//configuration file
		   break;
				//command line
	    case 1:
// This variant for console - April 2002
			if(strstr(main_name,"main_name.seq")!=NULL) 
			 strcpy(main_name,argv[i]);
			else flagname=1; //first name occupied yet

		   if((flagname)&&(strstr(second_name,"second_name.se2")!=NULL) )
			strcpy(second_name,argv[i]);

		   break;
	    default:fprintf(stdout,"Command line key not recognised. Please, check program parameters!\n");
	    break;
		}
	 }
   use_forward=(use_direct_search||use_complementary_search);
   use_reverse=(use_inverse_search||use_symmetry_search);
			
   i=0;while( (i<(int)strlen(main_name))&&(main_name[i]!='.') )
   {FASTAtable_name[i]=main_name[i];i++;}
    strcpy(FASTAtable_name+i,".out");
//char SignalDiscoveryFilename[50]="orlov_";
j=(int)strlen("orlov__")-1;
i=0; while( (i<(int)strlen(main_name))&&(i<50)) 
{SignalDiscoveryFilename[j+i]=main_name[i];i++;}
j=(int)strlen("orlov2__")-1;
i=0; while( (i<(int)strlen(main_name))&&(i<50)) 
{SignalDiscoveryFilename2[j+i]=main_name[i];i++;}

if(!(use_forward||use_reverse) ) 
{
use_direct_search=1;
use_forward=1;
fprintf(std_msg,"Repeat search not defined - Using default direct copying\n");
//return 0;
}

// fin=fopen(do_path(task_path,main_name),"rb");
 fin=fopen(main_name,"rb");
//if(count_method==2) 
// fsec=fopen(do_path(task_path,second_name),"rb");

     {
     buff0=new byte[SEQUENCE_LEN];
     if(buff0!=NULL)
       {
       base=new genfileinput();
       base->input_init(alphabet_pointer);
       if(base->success_init())
	 {

		 
		ALPHABET_LEN=base->get_alphabet_len();

		j=(int)(log((double)FRAME_LEN)/log((double)ALPHABET_LEN));//frame length
		//if(NSIZE<j+1) 
			NSIZE=j+1;//optimal
		if(NSIZE>FRAME_LEN)NSIZE=FRAME_LEN;
			
		i=(int)(log((double)(262144))/log((double)ALPHABET_LEN));//operating memory
		if(NSIZE>i) NSIZE=i;
		if(NSIZE<1) NSIZE=1;//minimal requirement

		if(debugregime==1)
		{
		//j=1; for(i=0;i<NSIZE;i++) j*=ALPHABET_LEN;//corrected for operating memory

        fprintf(std_msg,"ALPHABET_LEN=%i, j=logN(FRAME)=%i i=%i(memory<250K), set NSIZE=%i\n",ALPHABET_LEN,j,i,NSIZE);
        fprintf(std_msg,"count method=%i\n",count_method);
	fprintf(std_msg,"debugregime=%i\n",debugregime);
	fprintf(std_msg,"how=?? FRAME_LEN=%i SLIDE_STEP=%i\n",FRAME_LEN,SLIDE_STEP);
	}

		if(MarkovModelOrder>NSIZE) 
		{MarkovModelOrder=NSIZE;fprintf(std_msg,"Corrected MarkovModelOrder=%i\n",MarkovModelOrder);}
		//j=1; for(i=0;i<NSIZE;i++) j*=ALPHABET_LEN;
		j=(int)(log(8000.)/log((double)ALPHABET_LEN));//for entropy estimations
		if(MarkovModelOrder>j) {MarkovModelOrder=j;fprintf(std_msg,"Corrected MarkovModelOrder=%i\n",MarkovModelOrder);}

	 if((!best_function)&&(deffunc[0]))
	   {
	   if(!set_complementary_function())fprintf(std_msg,"Bad complementary function\n");
	   }
//	 if(best_function) work=new multidivizion(NSIZE,ALPHABET_LEN);
//	 else
//		{
//		if((count_method==0)&&(SLIDE_STEP==1))work=new slide_divizion(NSIZE,ALPHABET_LEN);
//		  else 
			  work=new divizion(NSIZE,ALPHABET_LEN);

//		}

	 if(work->success_init())
	   {
	   work->set_nonexact(0);//nonexact------ add here Feb.2004
	   if(alphabet_legend[0]) //something was defined by user
	   {
	   if(base->set_legend(alphabet_legend)) 
		   fprintf(std_msg,"bad user defined legend %s, use auto defined %s (from %s)\n",alphabet_legend,base->legend,alphabet_pointer);
	   }
	     if(seq_init())
	       {
	     switch(count_method) //main calculation of complexity profile
		   {
		   case 0:calculateFASTAprofile(fin,fout);
			  break;
		   case 1:calculateLCW(fin,fout);
			  break;
		   //case 1:calculate_whole(fin,fout);
			 // break;
//		   case 5:calculateFASTAprofile(fin,fout);
//			  break;
		   default:
			   fprintf(std_msg,"Method is not defined/not recognised from command line\n");
            break;
		 }//switch
	      }else fprintf(std_msg,"Not enought memory for decompose subsystem1\n");

		}else fprintf(std_msg,"Not enought memory for decompose subsystem2\n");

	 work->done();
	 delete work;
	 }
	   else fprintf(std_msg,"Not enought memory for file subsystem\n");

       delete base;
       delete buff0;
       }else fprintf(std_msg,"Not enought memory for buffer\n");
     if(fin!=NULL && fin!=stdin)fclose(fin);
     if(fout!=NULL && fout!=stdout)fclose(fout);
     }
   }else fprintf(std_msg,"No parameters in command line. In console: Type .exe -h for help\n");

 if(std_msg!=NULL && std_msg!=stdout)fclose(std_msg);
 return 0;
}

//Otstoy---------- reserve was below
