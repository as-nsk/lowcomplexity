#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ginput2.h"
#include "divizion2.h"
#include "tandem.h"
#define SEQUENCE_LEN 10000000 //100Mb
#define MAX_FRAME_LEN SEQUENCE_LEN
#define MIN_FRAME_LEN 5
#define DEFAULT_FRAME_LEN 50
#define MIN_SLIDE_STEP 1
#define MAX_N_SPEC_STAT 100 //100
#define MAX_N_SPEC_STAT2 1500//2000
#define MAX_SEQUENCIES 100//1024 #define MaxNumberofSeqs 200 //1024
//char default_alphabet[30]="atgc";
#define MAXDISTARRAYSIZE 10 //array of distanses by [x;x+DISTPARAMETER]
char protein_alphabet[21]="ACDEFGHIKLMNPQRSTVWY";
int FRAME_LEN=DEFAULT_FRAME_LEN;	// 50 step 1
int SLIDE_STEP=MIN_SLIDE_STEP;
int stdin_flag=0;
int show_help=0;
int detail_report=0;
int econom_report=0;
int special_report=0;
int tandem_report=0;
int user_max_level=20; //updated in April2002
int user_min_level=1;
int special_stat[MAX_N_SPEC_STAT][5]; //0-number, 1,2,3,4 - copy methods [0]=[1]+[2]+[3]+[4]
int special_stat_long[MAX_N_SPEC_STAT2+1][6]; // [][0=component_len,1=num_components,2,3,4,5 - methods]
int special_statl2[MAX_N_SPEC_STAT2+1];
int special_dist[MAXDISTARRAYSIZE+1][5];
int special_dist20[MAXDISTARRAYSIZE+1][5];
int special_max_n[5]={0,0,0,0,0};//maximal repeats of all types
int cur_max_n_spec_stat2=0;
int DISTPARAMETER=100; //interval for distance between repeats calculation

int count_method=0;
int best_function=0;
int use_direct_search=0;
int use_inverse_search=0;
int use_symmetry_search=0;
int use_complementary_search=0;
int use_reverse=0;
int use_forward=0;
int whole_sequencies=0;
int NSIZE=9;//9
int ALPHABET_LEN=4;
char alphabet_pointer[30]="atgc";//=default_alphabet;
static char main_name[64]="main_name.seq";
char second_name[64]="second_name.se2";
char res_name[64]="*.seq";
char rep_name[64]="*.rep";
char deffunc[MAX_ALPHABET_LEN+1]="tacg";
byte FR[MAX_ALPHABET_LEN]={1,0,3,2,4,5,6,7,8,9,10,11};//ATGC->TACG
byte FN[MAX_ALPHABET_LEN]={0,1,2,3,4,5,6,7,8,9,10,11};//ATGC->ATGC
char *abbrev[5]={"NW","D>","S<","C>","I<"};
char *hlp_text[]={"Decompose utility",
		  "Use .exe filename.ext [switches]",
		  "These switches available",
		  "<Switches>",
		  "   -aXX...: set alphabet (-aatgc default)",
		  "   -b: find best function by other step",
		  "   -c[disc]: set search&comparsion methods (-cdisc default)",
		  "      Letter d,i,s,c means these search methods:",
		  "      Direct, Inverse, Symmetry & Complementary search methods",
		  "   -zXX...: set func (default -dtacg, no effect with -b)",
		  "   -fXX set frame size (50 default), use only with -mp",
		  "   -tXX set slide step (50 default), use only with -mp",
		  "   -wXX for special report only - step for histogram of repeat distances (w=100 default)",
		  "   -h or -?: show this text",
		  "   -i: get data from standart input",
		  "//   -k(filename): machine configuration file",
		  "   -m(p|w|a|t|v|f): set count method",
		  "      These methods available:",
		  "      -mp count profile",
		  "      -mw decompose whole file",
		  "      -ma analyse file by other file: filename.se2",
		  "      -mt analyse all sequences by all others sequences",
		  "      -mv count complexity vector (can't use with -b)",
		  "      -mf count profiles for set in FASTA format",
		  "   -r(n|f|e|s) No report, Full detailed report, Economic or Special",
		  "   -y min complexity level for economic output",
		  "   -x max complexity level for economic output",
		  "//   -n(filename): set second file, use only with -ma",
		  "//   -o(filename): set output filename (default *.seq)",
		  "//   -r(filename): show detail information to filename",

		  "Example: .exe filename.seq",
		  "Attention: input data must be in FASTA format",
		  ""};
//******************************
int FASTAreport=1;
int MaxNumSeqs=MAX_SEQUENCIES;//was 1000
int MaxNumSteps=600;

 int TableDist[MAX_SEQUENCIES][MAX_SEQUENCIES];
 int TableDist2[MAX_SEQUENCIES][MAX_SEQUENCIES];
 double TableDistNormalised[MAX_SEQUENCIES][MAX_SEQUENCIES];
 double TableDistNormalised2[MAX_SEQUENCIES][MAX_SEQUENCIES];
 int TableDist3[MAX_SEQUENCIES][MAX_SEQUENCIES];
 double TableDist4[MAX_SEQUENCIES][MAX_SEQUENCIES];

//-- distances between repeats --
double special_dist_sum1[6];
double special_dist_qsum1[6];
double special_dist20_sum1[6];
double special_dist20_qsum1[6];
int special_dist_number=0;


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
	 case 'a':for(i=0;(i<100)&&(key[i+2]);i++)alphabet_pointer[i]=key[i+2];
		  alphabet_pointer[i]=0;
		  break;
	 case 'b':if(key[2]=='0')best_function=0;
		    else best_function=1;
		  break;
//	 case 'c':i=2; old - disc
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
	 case 'w':
		 DISTPARAMETER=atoi(key+2);
		  if(DISTPARAMETER<0) DISTPARAMETER=1;
		  break;
     case 'z':if(strlen(key+2)>2) strcpy(deffunc,key+2); //at least 2 symbols for complementary
				else fprintf(std_msg,"Bad complementary function defined in command line. Use default: %s\n",deffunc);
		 break;
	 case 'h':
	 case '?':i=0;
		  show_help=1;
		  while(hlp_text[i][0])fprintf(std_msg,"%s\n",hlp_text[i++]);
		  break;
//	 case 'i':stdin_flag=1;
//		  break;
//	 case 'k':
//		  h=fopen(key+2,"rt");
//		  if(h!=NULL)handle_file(h);
//		  fclose(h);
//		  break;
	 case 'm':switch(lowcase(key[2]))
			{
			case 'p':count_method=0;
				 break;
			case 'w':count_method=1;
//				 if(lowcase(key[3])=='s')whole_sequencies=1;
//				   else whole_sequencies=0;
				 whole_sequencies=0;
				 break;
			case 's':count_method=1;
				 whole_sequencies=1;
				 break;
			case 'a':count_method=2;
				 break;
			case 't':count_method=3;
				 break;
			case 'v':count_method=4;
				 break;
			case 'f':count_method=5;
				 break;
			default:  		 
				 fprintf(std_msg,"Bad calculation method defined in command line. Use default: %d\n",count_method);
				 break;
			}
		  break;
//	 case 'n':strcpy(second_name,key+2);
//		  break;
//	 case 'o':strcpy(res_name,key+2);
//		  break;
	 case 'r':
		 //if(key[2])
		 //   {
		 //   strcpy(rep_name,key+2);
		 //  detail_report=1;
		 //   }else detail_report=0;
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
			case 't':
					tandem_report=1;
				 break;
			case 'n':
			default: 
				detail_report=0;
				econom_report=0;
				special_report=0;
				tandem_report=0;
				 break;
			} 
		 break;
	 case 'o':strcpy(res_name,key+2);
		  break;
	 case 'y': //define min (y) and max (x) for 'econom' output
			{ 
			user_min_level=atoi(key+2);
			if(user_min_level<0)user_min_level=0;
			if(user_min_level>MAX_FRAME_LEN)user_min_level=MAX_FRAME_LEN-1;
			}
		  break;
	 case 'x': 
			{ 
			user_max_level=atoi(key+2);
			if(user_max_level>MAX_FRAME_LEN)user_max_level=MAX_FRAME_LEN;
			if(user_max_level<0)user_max_level=1;
			}
		  break;
	  default:return 2;
	 }
   return 0;
   }else return 1;
}

//******************************
/*
int handle_file(FILE *it)
{char *keyword[]={"alphabet","method","report","frame","slide_step","base",
		  "output","input","function","best","compare",""};
 char *curr_methods[]={"profile","whole_file","analyse","table","vector","whole_sequencies",""};
 char msg[1024],*hlp;
 char k[64];
 int i,j;
 while(!feof(it))
      {
      fgets(msg,1024,it);
      if((msg[0]=='/')&&(msg[1]=='/'))continue;
      sscanf(msg,"%s",k);
      i=0;
      while(strcmp(k,keyword[i])&&(keyword[i][0]))i++;
      if(keyword[i][0])
	{
	hlp=msg+strlen(k);
	sscanf(hlp,"%s",k+2);// rewrite second word to keyword string(argv) (meaning of parameter in the second column)
	k[0]=SWITCH_CHAR;    // key-symbol, here "-"
	switch(i)
	      {
	      case 0:
	      case 2:
	      case 3:
	      case 4:
	      case 10:
	      case 6:k[1]=keyword[i][0];//alphabet, frame, slide, compare
		     break;
	      case 7:strcpy(main_name,k+2);
		     break;
	      case 8:k[1]='d';//set function
		     break;
	      case 9:k[1]='b';
		     break;
	      case 5:k[1]='n';//set base file
		     break;
	      case 1:k[1]='m';
		     j=0;
		     while(strcmp(curr_methods[j],k+2)&&(curr_methods[j][0]))j++;
		     if(curr_methods[j][0])k[2]=curr_methods[j][0];
		       else k[2]=curr_methods[0][0];
		     if(j==5)k[3]='s';
		       else k[3]=0;
		     k[4]=0;
		     break;
	      }
	handle_key(k);
	}
      }
 return 0;
}
*/
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

int compare_spec(const void *arg1, const void *arg2 );
int compare_spec(const void *arg1, const void *arg2 )
   {
   return ((int *)(arg1)-(int *)(arg2));
   }

//--add April
int show_special_statistics(FILE *repfile);
int show_special_statistics(FILE *repfile)
{
char *spec_out_comment[6]={"Total   ","Direct  ","Symmetr.","Complem.","Inverted",""};
int i,j;//,j2sort;
int i1,i2,imin,tmp,minval;

double total_comp[5],total_len[5],mean_len[5];
double total_comp_more21[5],total_len_more21[5],mean_len_more21[5];
double aver_dist1[5];
double disp_dist1[5];
double long_total_comp[5], long_total_len[5], long_mean_len[5];

for(j=0;j<5;j++) 
	{
	total_comp[j]=0;total_len[j]=0;mean_len[j]=0;
	total_comp_more21[j]=0;total_len_more21[j]=0;mean_len_more21[j]=0;
	long_total_comp[j]=0;long_total_len[j]=0;long_mean_len[j]=0;
	aver_dist1[j]=0;disp_dist1[j]=0;
	}
//-define maximal repeat and its type
i1=special_max_n[1];i2=1;
for(j=1;j<5;j++) 
	if(special_max_n[j]>i1) 
	{i1=special_max_n[j];i2=j;}
fprintf(repfile,"Maximal Repeats = %i (type) %s:\n",i1,spec_out_comment[i2]);
for(j=1;j<5;j++) fprintf(repfile,"\t%i",special_max_n[j]);
fprintf(repfile,"\n");
fprintf(repfile,"Repeats ---- Type of repeats:\n");
fprintf(repfile,"Length Number    1-Standard 2-Symm. 3-Compl. 4-Inverted, (Check 1+2+3+4=Number)\n");

for(i=0;i<MAX_N_SPEC_STAT;i++)
{
fprintf(repfile,"%i\t%i\t",i,special_stat[i][0]);
for(j=1;j<5;j++) fprintf(repfile,"\t%i",special_stat[i][j]); //type of repeat
fprintf(repfile,"\n");

if(special_stat[i][0]>0)
	{
	for(j=0;j<5;j++) total_comp[j]+=special_stat[i][j];
	for(j=0;j<5;j++) total_len[j]+=((double)i)*((double)(special_stat[i][j]));
	}
}
//--------Excel log histogram August 2003
fprintf(repfile,"Excel log histogram August 2003\n");
/*for(i=1;i<MAX_N_SPEC_STAT;i++)
{
fprintf(repfile,"%i\t%i\t",i,special_stat[i][0]);
for(j=1;j<5;j++) fprintf(repfile,"\t%.2f",log((double)(special_stat[i][j])+0.5)/log(2.)); //type of repeat
fprintf(repfile,"\n");
}
fprintf(repfile,"--------------\n");
*/
for(i=1;i<MAX_N_SPEC_STAT;i++)
{
fprintf(repfile,"%i\t%i\t",i,special_stat[i][0]);
for(j=1;j<5;j++) fprintf(repfile,"\t%.2f",log((double)(special_stat[i][j])+1)/log(2.)); //type of repeat
fprintf(repfile,"\n");
}
fprintf(repfile,"--------------\n");


//------- statistics - more 21
for(i=22;i<MAX_N_SPEC_STAT;i++)
{
if(special_stat[i][0]>0)
	{
	for(j=0;j<5;j++) total_comp_more21[j]+=special_stat[i][j];
	for(j=0;j<5;j++) total_len_more21[j]+=((double)i)*((double)(special_stat[i][j]));
	}
}

 fprintf(repfile,"Longer components\n");
for(i=0;i<MAX_N_SPEC_STAT2;i++)
if(special_stat_long[i][1]>0)
	{
    fprintf(repfile,"%i\t%i\t",special_stat_long[i][0],special_stat_long[i][1]);
    for(j=1;j<5;j++) fprintf(repfile,"\t%i",special_stat_long[i][1+j]); //type of repeat
    fprintf(repfile,"\n");

	for(j=0;j<5;j++) total_comp[j]+=special_stat_long[i][1+j];//	total_comp+=special_stat_long[i][1];
	for(j=0;j<5;j++) total_len[j]+=((double)(special_stat_long[i][0]))*((double)(special_stat_long[i][1+j])); //lengths for all
//--only for long components
	for(j=0;j<5;j++) long_total_comp[j]+=special_stat_long[i][1+j];//	total_comp+=special_stat_long[i][1];
	for(j=0;j<5;j++) long_total_len[j]+=((double)(special_stat_long[i][0]))*((double)(special_stat_long[i][1+j])); //lengths for all
	}
/*
for(i=0;i<MAX_N_SPEC_STAT2;i++)
if(special_stat_long[i][1]>0)
 special_statl2[i]=special_stat_long[i][0];
   qsort( (void *)special_statl2, (size_t)MAX_N_SPEC_STAT2, sizeof(int), compare_spec );
 fprintf(repfile,"test sorted -------------Longer components\n");
for(i=0;i<MAX_N_SPEC_STAT2;i++)
{
   fprintf(repfile,"%i\n",special_statl2[i]);
}
*/
 fprintf(repfile,"2 -- sorted -------------Longer components\n");


for(i1=0;i1<MAX_N_SPEC_STAT2;i1++)
{
imin=i1;
minval=special_stat_long[i1][0];
for(i2=i1+1;i2<MAX_N_SPEC_STAT2;i2++)
{
if(minval>special_stat_long[i2][0])
	{
	minval=special_stat_long[i2][0];
	imin=i2;
	}
}
for(j=0;j<6;j++) 
	{
	tmp=special_stat_long[i1][j];
	special_stat_long[i1][j]=special_stat_long[imin][j];
	special_stat_long[imin][j]=tmp;
	}
}

for(i=0;i<MAX_N_SPEC_STAT2;i++)
{
if(special_stat_long[i][1]>0)
	{
    fprintf(repfile,"%i\t%i\t",special_stat_long[i][0],special_stat_long[i][1]);
    for(j=1;j<5;j++) fprintf(repfile,"\t%i",special_stat_long[i][1+j]); //type of repeat
    fprintf(repfile,"\n");
	}
}


//Last unsorted (no place in array
 fprintf(repfile,"\nUnsorted list: length\t%i\t\n",special_stat_long[MAX_N_SPEC_STAT2][1]);
    for(j=1;j<5;j++) fprintf(repfile,"\t%i",special_stat_long[MAX_N_SPEC_STAT2][1+j]); //type of repeat
    fprintf(repfile,"\n");

// substitute  MAX_N_SPEC_STAT2 by maximal
 if( special_stat_long[MAX_N_SPEC_STAT2][1] > 0)
 {
 for(j=0;j<5;j++) total_comp[j]+=special_stat_long[MAX_N_SPEC_STAT2][1+j];
 for(j=0;j<5;j++) total_len[j]+=(double)(MAX_N_SPEC_STAT+MAX_N_SPEC_STAT2)*((double)(special_stat_long[MAX_N_SPEC_STAT2][1+j]));
 }

 for(j=0;j<5;j++) 
  if(total_comp[j]>0)
   mean_len[j]=total_len[j]/total_comp[j];
 for(j=0;j<5;j++) 
  if(long_total_comp[j]>0)
   long_mean_len[j]=long_total_len[j]/long_total_comp[j];
 
 for(j=0;j<5;j++) 
 {
 total_len_more21[j]+=long_total_len[j];
 total_comp_more21[j]+=long_total_comp[j];
 }
 for(j=0;j<5;j++) 
  if(total_comp_more21[j]>0)
   mean_len_more21[j]=total_len_more21[j]/total_comp_more21[j];

 fprintf(repfile,"Complexity = number of components %.2f\n",total_comp[0]);
 fprintf(repfile,"Length of sequence analysed - Check length=(sum of components) %.2f\n",total_len[0]);
//------------- just for paper - August 2003
FILE *comparetablefile;
 comparetablefile=fopen("table1.txt","a+");
 fprintf(comparetablefile,"%.2f\t%.2f\t%.2f\t%.2f\n",total_len[0],mean_len[0],log((double)(total_len[0])),log((double)(total_len[0]))/mean_len[0]);
fclose(comparetablefile);
//---------------
 fprintf(repfile,"Mean component length %.2f\n",mean_len[0]);
 fprintf(repfile,"Logarithm length %.2f\n",log((double)(total_len[0])) );
 fprintf(repfile,"Relation constant log (N) / Mean  %.2f\n",log((double)(total_len[0]))/mean_len[0] );
 fprintf(repfile,"Num. of component  Length(nt)  Mean\n");
//   1-Standard 2-Symmetry 3-Compl 4-Inverted, (Check 1+2+3+4=Number)
for(j=0;j<5;j++)
 {
 fprintf(repfile,"%s",spec_out_comment[j]);
 fprintf(repfile,"\t%i",(int)total_comp[j]);
 fprintf(repfile,"\t%i",(int)total_len[j]);
 fprintf(repfile,"\t%.3f\n",mean_len[j]);
 }

fprintf(repfile,"Longer (>100) components\n Number Length(nt) Mean\n");
 for(j=0;j<5;j++)
 {
 fprintf(repfile,"%s",spec_out_comment[j]);
 fprintf(repfile,"\t%.0f",long_total_comp[j]);
 fprintf(repfile,"\t%.0f",long_total_len[j]);
 fprintf(repfile,"\t%.3f\n",long_mean_len[j]);
 }

//--debug check1
 fprintf(repfile,"----- Check1 total length ----\n");
 fprintf(repfile,"special_dist_number=%i\n",special_dist_number);
 for(j=0;j<5;j++)
  fprintf(repfile,"total_comp[%i]=%.i\n",j,(int)total_comp[j]);

 fprintf(repfile,"\n Ratio of components (Longer>100)/Total: Number, Length /TotalLength\n");
 for(j=0;j<5;j++)
 {
 fprintf(repfile,"%s",spec_out_comment[j]);
 fprintf(repfile,"\t%.4f",long_total_comp[j]/total_comp[j]);
 fprintf(repfile,"\t%.4f",long_total_len[j]/total_len[j]);
 fprintf(repfile,"\t%.4f\n",long_total_len[j]/total_len[0]);
 }

fprintf(repfile,"\nLonger (>21) components\n Number Length(nt) Mean  >(21+100)/TotalLength\n");
 for(j=0;j<5;j++)
 {
 fprintf(repfile,"%s",spec_out_comment[j]);
 fprintf(repfile,"\t%.2f",total_comp_more21[j]);
 fprintf(repfile,"\t%.2f",total_len_more21[j]);
 fprintf(repfile,"\t%.2f",mean_len_more21[j]);
 fprintf(repfile,"\t%.2f\n",(total_len_more21[j]+long_total_len[j])/total_len[0]);
 }
//--debug check2
 fprintf(repfile,"----- Check2 total length ----\n");
 fprintf(repfile,"special_dist_number=%i\n",special_dist_number);
 for(j=0;j<5;j++)
  fprintf(repfile,"total_comp[%i]=%i\n",j,(int)total_comp[j]);

 fprintf(repfile,"---------- Distance D. S. C. I. -------------\n");
 for(i=0;i<=MAXDISTARRAYSIZE;i++)
 {
 fprintf(repfile,"%i",i*DISTPARAMETER);
 for(j=0;j<5;j++)
	 fprintf(repfile,"\t%i",special_dist[i][j]);
 fprintf(repfile,"\n");
 }
 fprintf(repfile,"------- Distance (Fragment>20) D. S. C. I. -----------\n");
 for(i=0;i<=MAXDISTARRAYSIZE;i++)
 {
 fprintf(repfile,"%i",i*DISTPARAMETER);
 for(j=0;j<5;j++)
	 fprintf(repfile,"\t%i",special_dist20[i][j]);
 fprintf(repfile,"\n");
 }

 fprintf(repfile,"----- Check total length ----\n");
 fprintf(repfile,"special_dist_number=%i\n",special_dist_number);
 for(j=0;j<5;j++)
  fprintf(repfile,"total_comp[%i]=%i\n",j,(int)total_comp[j]);

 fprintf(repfile,"----- Average distance ----\n");
 for(j=0;j<5;j++)
 {
 if(total_comp[j]>1)
	{
	aver_dist1[j]=special_dist_sum1[j]/(double)(total_comp[j]);
	disp_dist1[j]=(special_dist_qsum1[j]-(special_dist_sum1[j]*special_dist_sum1[j]/(double)(total_comp[j])))/(double)(total_comp[j]-1);
	disp_dist1[j]=sqrt(disp_dist1[j]);
	}
 }
 fprintf(repfile,"Aver:");
 for(j=0;j<5;j++)
	 fprintf(repfile,"\t%i",(int)aver_dist1[j]);
 fprintf(repfile,"\n");
 fprintf(repfile,"Disp:");
 for(j=0;j<5;j++)
	 fprintf(repfile,"\t%.2f",disp_dist1[j]);
 fprintf(repfile,"\n");

 fprintf(repfile,"----- Average distance (>21)----\n");
 for(j=0;j<5;j++)
 {
	aver_dist1[j]=0;
	disp_dist1[j]=0;
 if(total_comp[j]>1)
	{
	aver_dist1[j]=special_dist20_sum1[j]/(double)(total_comp_more21[j]+long_total_comp[j]);
	disp_dist1[j]=(special_dist20_qsum1[j]-((special_dist20_sum1[j]*special_dist20_sum1[j])/(double)(total_comp_more21[j]+long_total_comp[j]) ))/(double)((total_comp_more21[j]+long_total_comp[j])-1);
	disp_dist1[j]=sqrt(disp_dist1[j]);
	}
 }
 fprintf(repfile,"Aver:");
 for(j=0;j<5;j++)
	 fprintf(repfile,"\t%i",(int)aver_dist1[j]);
 fprintf(repfile,"\n");
 fprintf(repfile,"Disp:");
 for(j=0;j<5;j++)
	 fprintf(repfile,"\t%.2f",disp_dist1[j]);
 fprintf(repfile,"\n");


return(0);
}


int special_statistics(int n, int method, int ndist);
int special_statistics(int n, int method, int ndist)
{
int i,flag=0;
int kdist;
double tmpdouble;

if(n<0) return(-1); //error
if(ndist<0) return(-1); //error
if(method<0) return(-1); //error
if(method==0)method=1; //mark as direct repeat for NewLetter generation by default
for(i=1;i<=4;i++)if(n>special_max_n[method])special_max_n[method]=n;
if(n<MAX_N_SPEC_STAT) 
	{
	special_stat[n][0]++;
	special_stat[n][method]++; //methods=1,2,3,4 [n][0]=[n][1]+[n][2]+[n][3]+[n][4]
	}
else
{
i=0;
while ((i<cur_max_n_spec_stat2) && (!flag) )
	{
	if(special_stat_long[i][0]==n) 
	{special_stat_long[i][1]++; //total number of repeats
	special_stat_long[i][1+method]++;//number of repeat of this type
	flag=1;}
	i++;
	}

	if(!flag) 
	{
	if (cur_max_n_spec_stat2<MAX_N_SPEC_STAT2-1 ) //new record in this array
		{
		cur_max_n_spec_stat2++;
		special_stat_long[cur_max_n_spec_stat2][0]=n;			//fix new repeat length
		special_stat_long[cur_max_n_spec_stat2][1]=1;			//fix number of repeat  - now 1 only
		special_stat_long[cur_max_n_spec_stat2][1+method]=1;	//fix type of repeat
		}
	else 
		{
		special_stat_long[MAX_N_SPEC_STAT2][1]++;// no place in array, fix additional number of repeat
		special_stat_long[MAX_N_SPEC_STAT2][1+method]++;// fix type of repeat
		}
	}
}

//July 2003
    kdist=(int)( (double)(ndist)/(double)(DISTPARAMETER) );
    if(kdist>MAXDISTARRAYSIZE) kdist=MAXDISTARRAYSIZE;
//all repeats
	special_dist[kdist][0]++;
	special_dist[kdist][method]++;

tmpdouble=(double)ndist;//to calculate sum, dispersion of distances beetwen repeart
//-------- sum, dispersion
special_dist_sum1[0]+=tmpdouble;
special_dist_sum1[method]+=tmpdouble;
special_dist_qsum1[0]+=tmpdouble*tmpdouble;
special_dist_qsum1[method]+=tmpdouble*tmpdouble;

if(n>21)//longer repeats
	{
    special_dist20[kdist][0]++;
    special_dist20[kdist][method]++;
//-------- sum, dispersion
	special_dist20_sum1[0]+=tmpdouble;
	special_dist20_sum1[method]+=tmpdouble;
	special_dist20_qsum1[0]+=tmpdouble*tmpdouble;
	special_dist20_qsum1[method]+=tmpdouble*tmpdouble;
	}

special_dist_number++;
return(0);
}


int calculate_profile(FILE *f0,FILE *rep) //main profile calculation
{int result,how,from,size,number=0,complexity,i=0,j=0;
 list_info hdr={{0,0,0,0},NULL};
 list_info* resu=&hdr;
 info h0;
int FASTAcurrstep=0;//further - from 0 to number of sequences

// if(SLIDE_STEP>FRAME_LEN)SLIDE_STEP=FRAME_LEN;
	   if(SLIDE_STEP<MIN_SLIDE_STEP)SLIDE_STEP=MIN_SLIDE_STEP;
 if(detail_report)fprintf(std_msg,"Frame length:%i\nSlide step:%i\n",FRAME_LEN,SLIDE_STEP);

//temp/ August 2003
 int previouscomplexity;
 previouscomplexity=-FRAME_LEN;
 FILE *tandemfile;
tandemfile=fopen("economrep.txt","w");
 if(econom_report)
	{
	fprintf(std_msg,"Frame length:%i\nSlide step:%i\n",FRAME_LEN,SLIDE_STEP);
	fprintf(std_msg,"User-defined complexity less than Min or greater than Max\n Min:%i\n Max:%i\n",user_min_level,user_max_level);

	fprintf(tandemfile,"Frame length:%i\nSlide step:%i\n",FRAME_LEN,SLIDE_STEP);
	fprintf(tandemfile,"User-defined complexity less than Min or greater than Max\n Min:%i\n Max:%i\n",user_min_level,user_max_level);
	}

 base->set_file(f0);
//--------------- main cycle
 do{
   result=fill_sequence(buff0,&how); //find and set next sequence in the file (in FASTA format)
   if(!how)continue;

   if(detail_report)
     {
     fprintf(std_msg,"Sequence No %i\n",number+1);
     clean_statistics();
     }
   if(econom_report) 
	{
     fprintf(std_msg,"Sequence No %i\n",number+1);
     fprintf(tandemfile,"Sequence No %i\n",number+1);
     clean_statistics();
	}
/*
   if(FASTAreport) //sequence number = "number", current step (position) = 0
	{
	FASTAcurrstep=0;
	}
*/
   from=0;
   if((!best_function)&&(SLIDE_STEP==1))
     {
     i=j=0;
     resu=&hdr;
     ft_vector2=0;
     complexity=0;
     work->set_universe(buff0,how);
     if(how<FRAME_LEN)size=how;
       else size=FRAME_LEN;

	while(i<size) //decompose first window in the sequence!
	{
	if(j<=i) //make component list for output
		{
		resu->next=new list_info;
		resu=resu->next;
		resu->dt=find_maximal_match(buff0+i,size-i,i,i);
		complexity++;
		resu->next=NULL;
		if(!resu->dt.method)resu->dt.from=i;
		resu->dt.seq_number=number;
		j+=resu->dt.len;
		}
	work->add_vector(buff0+i,i,size-i);
	i++;
	}
     from=1;
     fprintf(rep,"%i\t%i\n",from,complexity);//only profile without details
     collect_statistics(complexity);
/*   if(FASTAreport) 
	{
	FASTAtable[number][FASTAcurrstep]=complexity;
    if(FASTAcurrstep<MaxNumSteps) FASTAcurrstep++;
	//    FASTAcollect(number,currstep,complexity);
	}
*/
     if(detail_report)
	   {
	   show_item_list(hdr.next,buff0,0);//show list of components in window
	   }
    if(econom_report) 
		{
		if(complexity<user_min_level || complexity>user_max_level)
if(from-previouscomplexity>FRAME_LEN)
{
	previouscomplexity=from;
	   	   {
          fprintf(tandemfile,"low complexity %i seq. from %i by length %i\n",complexity, from,FRAME_LEN);
	      for(i=from;i<from+FRAME_LEN+FRAME_LEN;i++)fprintf(tandemfile,"%c",base->legend[buff0[i]]);
          fprintf(tandemfile,"\n");
			show_item_list(hdr.next,buff0,0);
			}
}
		}

	 while((from+FRAME_LEN)<how)//<=  decompose second etc. windows in the sequence!
	  {
//printf("number=%i, from=%i,ft_vector=%i\n",number,from,ft_vector2);
		 ((slide_divizion*)work)->slide_tree(FRAME_LEN);
	  hdr.next=new_decomposition(hdr.next,&complexity,number);
     fprintf(rep,"%i\t%i\n",from,complexity);//only profile without details?
		collect_statistics(complexity);
/*   if(FASTAreport) 
	{
	FASTAtable[number][FASTAcurrstep]=complexity;
    if(FASTAcurrstep<MaxNumSteps)	FASTAcurrstep++;
	}
*/
	  if(detail_report)
	    {
	    show_item_list(hdr.next,buff0,from);
	    }
     if(econom_report) 
		{
		if(complexity<user_min_level || complexity>user_max_level)
	   	   {
if(from-previouscomplexity>FRAME_LEN)
{
	previouscomplexity=from;
          fprintf(tandemfile,"low complexity seq. from %i by length %i\n",from,FRAME_LEN);
	      for(i=from;i<from+FRAME_LEN+FRAME_LEN;i++)fprintf(tandemfile,"%c",base->legend[buff0[i]]);
          fprintf(tandemfile,"\n");

		  show_item_list(hdr.next,buff0,0);
}
			}
	    }
		from++;
	  }
	 work->clear();				//clear tree for the next sequence
	 list_info *hlp,*wo=hdr.next;
	 while(wo!=NULL)
		  {
		  hlp=wo->next;
		  delete wo;
		  wo=hlp;
		  }
	 if(detail_report)show_statistics(std_msg);
	 if(econom_report)show_statistics(std_msg);
     }else						//Fillipov symbol iterative procedure if step!=1
     {
     do{
       if((from+FRAME_LEN)<=how)size=FRAME_LEN;
	 else size=how-from;
       work->clear();
       work->set_universe(buff0+from,size);
       int limit=0;
       complexity=0;
       i=0;

       do{
	 if(i>=limit)
	   {
	   h0=find_maximal_match(buff0+from+i,size-i,limit,limit);
	   h0.seq_number=number;

/*	   if(detail_report)
	     {
	     h0.from+=from;
	     show_item(h0,buff0+from+i,i+from);
	     }
     if(econom_report) 
		{
		if(complexity<user_min_level || complexity>user_max_level)
	   	   {
          fprintf(tandemfile,"low complexity seq. from %i by length %i\n",from,FRAME_LEN);
	      for(i=from;i<from+FRAME_LEN;i++)fprintf(tandemfile,"%c",base->legend[buff0[i]]);
          fprintf(tandemfile,"\n");
			show_item_list(hdr.next,buff0,0);
			}
		}
*/	   limit+=h0.len;
	   complexity++;
	   }
	 work->add_vector(buff0+from+i,i,size-i);
	 i++;
	 }while(i<size);

       from+=SLIDE_STEP;
       if(detail_report)
	     {
	     fprintf(std_msg,"Complexity is %i\n",complexity);
	     }
       if(econom_report) //here the same output as for detail
	     {
//	     fprintf(std_msg,"Complexity is %i\n",complexity);
//	     fprintf(tandemfile,"Complexity is %i\n",complexity);
		if(complexity<user_min_level || complexity>user_max_level)
	   	   {
if(from-previouscomplexity>FRAME_LEN)
{
	previouscomplexity=from;
          fprintf(tandemfile,"low complexity %i seq. from %i by length %i\n",complexity, from,FRAME_LEN);
	      for(i=from;i<from+FRAME_LEN+FRAME_LEN;i++)fprintf(tandemfile,"%c",base->legend[buff0[i]]);
          fprintf(tandemfile,"\n");

		  show_item_list(hdr.next,buff0,0);
}
			}
	     }
     fprintf(rep,"%i\t%i\n",from,complexity);//only profile without details
		 collect_statistics(complexity);
       }
	   while((from+size)<how);//here was < 

	 if(detail_report)show_statistics(std_msg);
	 if(econom_report)show_statistics(std_msg);
     }//end else
   number++;
   }while(result!=END_OF_DATA);

   fclose(tandemfile);//temp. long component output
 return 0;
}


//------------------------------------------------------------------------
int calculate_whole(FILE *f0,FILE *rep)
{int result,hows=0,i,j,limit,number=0;
int tempshow=0;
 info h0;
 //rep=rep;
//-----add for long component output
 FILE *tandemfile;
tandemfile=fopen("CurrentSpecialReport.txt","w");
//int i2,j2;

 if(detail_report)
   {
   if(whole_sequencies)fprintf(std_msg,"Compose whole sequencies:%s\n",main_name);
     else fprintf(std_msg,"Compose whole file:%s\n",main_name);
   }
 if(econom_report)
   {
	fprintf(std_msg,"User-defined complexity less than Min or greater than Max\n Min:%i\n Max:%i\n",user_min_level,user_max_level);
   if(whole_sequencies)fprintf(std_msg,"User-limited statistics output. Compose whole sequencies:%s\n",main_name);
   else 
 	fprintf(std_msg,"User-limited ouput. Compose whole file:%s\n",main_name);
   }

 base->set_file(f0);
 do{
   result=fill_sequence(buff0+hows,&i);
   if(!i)continue;
   add_sequence(hows);

  if(whole_sequencies)clean_statistics();
   fprintf(std_msg,"Sequence No %i (Length %i)\n",number+1,(int)i);

   work->set_universe(buff0,hows+i);
   limit=0;
   j=0;
   do{
     if(j>=limit)
       {
       h0=find_maximal_match(buff0+hows+j,i-j,limit+hows,limit+hows);
       if(detail_report)
	     {
		 if(whole_sequencies)h0.seq_number=number;
	       else h0=determine_sequence(h0);
		 collect_statistics(h0.len);
		 show_item(h0,buff0+hows+j,j);//(h0,buff0+hows+j,j+hows)
	     }
       if(special_report)
	     {
         special_statistics(h0.len,h0.method,(j-h0.from));//,h0.from,h0.method,j);//analyse repeat types j-current position?
	     }

//----Aug 2003here one maximal component output
       if(econom_report)
		{
		if(h0.len>user_max_level)
			{
			show_long_repeats(tandemfile, buff0, h0, j, limit, hows);
/*
//		 if(whole_sequencies)h0.seq_number=number;
//	       else h0=determine_sequence(h0);
          fprintf(tandemfile,"length=%i, pos=%i (lim=%i,from=%i,diff=%i), met=%i\n",h0.len,j, limit, h0.from,(j-h0.from),h0.method);
		  show_item(h0,buff0+hows+j,j);//(h0,buff0+hows+j,j+hows)
		  for(j2=0;j2<h0.len;j2++)
			  fprintf(tandemfile,"%c",base->legend[buff0[hows+j+j2]]);
          fprintf(tandemfile,"\n");
		i2=findminrepunit(buff0+hows+j,h0.len);
          if(i2<h0.len)fprintf(tandemfile,"Inside long fragment->Tandem - minimal rep.unit (%i bp *%.1f):",i2,(double)(h0.len)/(double)(i2));
		  else fprintf(tandemfile,"Long fragment - no simple repeat (%i bp=%i): ",i2,(int)h0.len);
		  for(j2=0;j2<i2;j2++)
	          fprintf(tandemfile,"%c",base->legend[buff0[hows+j+j2]]);
          fprintf(tandemfile,"\n");
//		  for(i2=0;i2<h0.len;i2++)
//			  fprintf(tandemfile,"%c",base->legend[buff0[hows+j+i2]]);
//          fprintf(tandemfile,"\n");
			*/
			}
	   }

	   if(tandem_report)
	   {
show_tandem_report(tandemfile, buff0, h0, j, limit, hows);//show info
/*
		if((limit<=h0.from-1) && ((h0.method==2)||(h0.method==4)) ) //find reversed repeat
		{
        fprintf(tandemfile,"Palindrome length=%i, pos=%i (lim=%i,from=%i,diff=%i), met=%i\n",h0.len,j, limit+hows, h0.from,(j-h0.from),h0.method);
		for(j2=0;j2<(limit-h0.from)+2*h0.len;j2++)
		  fprintf(tandemfile,"%c",base->legend[buff0[h0.from-h0.len+j2]]);
          fprintf(tandemfile,"\n");
		for(j2=0;j2<h0.len;j2++)
		  fprintf(tandemfile,"<");
		for(j2=0;j2<(limit-h0.from);j2++)
		  fprintf(tandemfile,"-");
		for(j2=0;j2<h0.len;j2++)
		  fprintf(tandemfile,">");
          fprintf(tandemfile,"\n");

		  i2=findminrepunit(buff0+hows+j,h0.len);
          if(i2<h0.len)fprintf(tandemfile,"Reverse->Tandem - minimal rep.unit (%i bp *%.1f):",i2,(double)(h0.len)/(double)(i2));
		  else fprintf(tandemfile,"Reverse palindrome, no rep.unit (%i bp=%i): (",i2,(int)h0.len);
		  for(j2=0;j2<i2;j2++)
	          fprintf(tandemfile,"%c",base->legend[buff0[hows+j+j2]]);
          fprintf(tandemfile,")\n");

		}

		if((h0.len>(limit-h0.from))&&(h0.len>1)&& ((h0.method==1)||(h0.method==3))) //find forward (tandem) repeat
{          
show_item(h0,buff0+hows+j,j+hows);//(h0,buff0+hows+j,j+hows)
          fprintf(tandemfile,"Overlap length=%i, pos=%i (lim=%i,from=%i,diff=%i), met=%i\n",h0.len,j, limit, h0.from,(j-h0.from),h0.method);
//fprintf(tandemfile,"---------------\n");
	if((hows+j-(limit-h0.from+5))>0)
			for(j2=0;j2<(limit-h0.from+5);j2++)
			  fprintf(tandemfile,"%c",base->legend[buff0[hows+j-(limit-h0.from+5)+j2]]);
          //fprintf(tandemfile,"[");
		  for(j2=0;j2<h0.len;j2++)
			  fprintf(tandemfile,"%c",base->legend[buff0[hows+j+j2]]);
          //fprintf(tandemfile,"]");
		  for(j2=0;j2<5;j2++)
			  fprintf(tandemfile,"%c",base->legend[buff0[hows+j+h0.len+j2]]);
          fprintf(tandemfile,"\n");
		  //mark coincindence
	if((hows+j-(limit-h0.from+5))>0)
			  for(j2=0;j2<(limit-h0.from+5);j2++) fprintf(tandemfile,"|");
		  for(j2=0;j2<h0.len;j2++)  fprintf(tandemfile,">");
		  for(j2=0;j2<5;j2++)  fprintf(tandemfile,"|");
          fprintf(tandemfile,"\n");
	if((hows+j-(limit-h0.from+5))>0)
		  for(j2=0;j2<5;j2++)
			  fprintf(tandemfile,"%c",base->legend[buff0[h0.from-5+j2]]);
		  for(j2=0;j2<h0.len;j2++)
	          fprintf(tandemfile,"%c",base->legend[buff0[h0.from+j2]]);
	if((hows+j-(limit-h0.from+5))>0)
		  for(j2=0;j2<5;j2++)
	          fprintf(tandemfile,"%c",base->legend[buff0[h0.from+h0.len+j2]]);
          fprintf(tandemfile,"\n");
	if((hows+j-(limit-h0.from+5))>0)
		  for(j2=0;j2<5;j2++) fprintf(tandemfile,"-");
		  for(j2=0;j2<h0.len;j2++)  fprintf(tandemfile,"*");
		  for(j2=0;j2<5;j2++)  fprintf(tandemfile,"-");
          fprintf(tandemfile,"\n");

		i2=findminrepunit(buff0+h0.from,h0.len);
          if(i2<h0.len)fprintf(tandemfile,"Tandem - minimal rep.unit (%i bp *%.1f):",i2,(double)(h0.len)/(double)(i2));
		  else fprintf(tandemfile,"Artefact, rep.unit (%i bp=%i): (",i2,(int)h0.len);
		  for(j2=0;j2<i2;j2++)
	          fprintf(tandemfile,"%c",base->legend[buff0[h0.from+j2]]);
          fprintf(tandemfile,")\n");
}*/
	   }

	   limit+=h0.len;
       // if(!special_report) 
//	   if(limit>10000*tempshow)  {fprintf(rep,"%i<-%i <%i>m=%i\n",limit,h0.from,h0.len,h0.method);tempshow++;}
       }
     work->add_vector(buff0+hows+j,hows+j,i-j);
     j++;
     }while(j<i);
   buff0[hows+i]=NEW_SEQUENCE_ITEM;
   hows+=(i+1);
   number++;

//   printf("Length(limit,i,j)=%i,%i,%i\n",limit,i,j);

 if(whole_sequencies)
     {
     work->clear();
     hows=0;
     }

   if(detail_report)show_statistics(std_msg);
   }while(result!=END_OF_DATA && hows<SEQUENCE_LEN);

 if(special_report) 
	{
	 show_special_statistics(rep);
	fprintf(rep,"Length(limit,i,j)=%i,%i,%i\n",limit,i,j);
	}
 if(tandem_report) 
	{
	fprintf(tandemfile," Summary tandem report: \n");
	fprintf(tandemfile," Reverted: Tandem= %i \tPalindrome= %i\n",revtandem,revpalindrome);
	fprintf(tandemfile," Direct:   Tandem= %i \tDir. close= %i\n",dirtandem,dirclose);
	}
 if(econom_report) 
	{
	fprintf(tandemfile," Summary long repeas: \n");
	fprintf(tandemfile," LongSSR = %i (%.2f)\n",LongSSR,(double)LongSSR/(double)(LongSSR+LongUniq));
	fprintf(tandemfile," LongUniq= %i (%.2f)\n",LongUniq,(double)LongUniq/(double)(LongSSR+LongUniq));
	}
 work->clear();

fclose(tandemfile);
 return 0;
}


//------------- 01.08.2003 - addition -no output, only sequence complexity -----------------------------------------------------------
int calc_whole_seqs(FILE *f0,FILE *rep)
{int result,hows=0,i,j,limit,number=0;
int tempshow=0;
int icomplexity;
 info h0;
 rep=rep;

 base->set_file(f0);
 do{
   result=fill_sequence(buff0+hows,&i);
   if(!i)continue;
   add_sequence(hows);
//   clean_statistics();
   work->set_universe(buff0,hows+i);
   icomplexity=0; //compexity for each sequence in file
   limit=0;
   j=0;
   do{
     if(j>=limit)
       {
       h0=find_maximal_match(buff0+hows+j,i-j,limit+hows,limit+hows);
	   limit+=h0.len;
	   icomplexity++;
		}
     work->add_vector(buff0+hows+j,hows+j,i-j);
     j++;
     }while(j<i);

   buff0[hows+i]=NEW_SEQUENCE_ITEM;
   hows+=(i+1);
   TableDist[number][number]=icomplexity;
   number++;

   fprintf(rep,"Seq.%i Complexity=%i Length(limit,i,j)=%i,%i,%i\n",number, icomplexity, limit,i,j);

// if(whole_sequencies)
     {
     work->clear();
//     hows=0;
     }
   }while(result!=END_OF_DATA);

   fprintf(rep,"Complexity of every sequence by itself is calculated %\n",number);

 return 0;
}


int calculate_analyse(FILE *f0,FILE *f1,FILE *rep)
{int result,i,j,how,hows=0,number=0;
 info h0;
 base->set_file(f0);
 if(detail_report)fprintf(std_msg,"Composed file:%s\nBase file:%s\n",main_name,second_name);
 do{
   result=fill_sequence(buff0+hows,&i);
   if(!i)continue;
   add_sequence(hows);
   for(j=0;j<i;j++)work->add_vector(buff0+hows+j,hows+j,i-j);
   if(result==NEW_SEQUENCE)
     {
     buff0[i+hows]=NEW_SEQUENCE_ITEM;
     i++;
     }
   hows+=i;
   }while(result!=END_OF_DATA);
 base->set_file(f1);
 i=0;
 do{
   result=fill_sequence(buff1,&how);
   i=0;
   work->set_universe(buff0,hows);
   while(i<how)
	{
	h0=find_maximal_match(buff1+i,how-i,hows,hows);
	if(detail_report)
	  {
	  h0=determine_sequence(h0);
	  show_item(h0,buff1+i,i);
	  }
	i+=h0.len;
	fprintf(rep,"%i\n",h0.len);
	}
   if(detail_report)fprintf(std_msg,"Sequence No %i\n",number+1);
   number++;
   }while(result!=END_OF_DATA);
 return 0;
}
//April 2003

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

int calculate_table(FILE *f0,FILE *rep)//each sequence by each sequence in the sample
{int *seq,i,j,howseq=0,how=0,result,complexity;
 int seqlengths[MAX_SEQUENCIES];
char astr[1024];
int limit,t,size,size1;
info h0;
int porog=20,com2num,com2len;

 seq=new int[MAX_SEQUENCIES];// was 1024

 Curr_seq_num=0;
 while(!feof(f0))
 {
 fscanf(f0,"%s",astr);
 if(astr[0]=='>')
 {
	 printf("\tOK %i\n",Curr_seq_num);
 for(j=0;j<20,j<(int)strlen(astr);j++) printf("%c",astr[j]);
 for(j=0;j<20,j<(int)strlen(astr);j++)Seq_names[Curr_seq_num][j]=astr[j+1];
 Curr_seq_num++;
 }
 }
 fseek(f0,0,SEEK_SET);
   fprintf(std_msg,"Scannned Curr_seq_num=%i\n",Curr_seq_num);


     fprintf(rep,"- TableDist - check -------------------\n");
 for(i=0;i<Curr_seq_num;i++)
	{
	fprintf(rep,"Seq%i\t",i);
	for(j=0;j<Curr_seq_num;j++)
		fprintf(rep,"%i ",TableDist[i][j]);
	fprintf(rep,"\n");
 }


 if(seq!=NULL)
   {
   for(i=0;i<MAX_SEQUENCIES;i++)seq[i]=0;
   base->set_file(f0);
   do{
     result=fill_sequence(buff0+how,&i);
     if(!i)continue;
     add_sequence(how);

	work->set_universe(buff0,how+i);
   complexity=0; //compexity for each sequence in file
   limit=0;
   com2num=0;
   com2len=0;
   j=0;
   do{
     if(j>=limit)
       {
       h0=find_maximal_match(buff0+how+j,i-j,limit+how,limit+how);
	   limit+=h0.len;
		  complexity++;
          if(h0.len>porog) com2num++;
          if(h0.len>porog) com2len+=h0.len;
		}
     work->add_vector(buff0+how+j,how+j,i-j);
     j++;
     }while(j<i);
	     TableDist[howseq][howseq]=complexity;//j-th decomposed by i-th
	     TableDist2[howseq][howseq]=com2num;
	     TableDist4[howseq][howseq]=(double)com2len/(double)(i);

     seq[howseq++]=how;
     how+=i;
	 	work->clear();
   }while(result!=END_OF_DATA && howseq<MAX_SEQUENCIES);

   seq[howseq]=how;
   if(detail_report)fprintf(std_msg,"Make table of decompose. %i sequences found.\n",howseq);

   for(i=0;i<howseq;i++)
   fprintf(std_msg,"%i\t%i(%i)\n",i,seq[i],seq[i+1]-seq[i]);


   if(howseq)
     {
     for(i=0;i<howseq;i++)
	{
	size=seq[i+1]-seq[i];
   seqlengths[i]=size;//add for table analysis - April

	if(detail_report)fprintf(std_msg,"Set %i sequence as a basic. Size is %i\n",i,size);
	fprintf(std_msg,"Seq%i[%i]\t",i+1,size);
//--------------- Self.calc.	
/*	work->clear();
	work->set_universe(buff0+seq[i],size);
	     complexity=0;
		 com2num=0;
		 com2len=0;
		 t=0;
		 j=0;
		 do
		 {
		  if(j>=t)
		  {
		  h0=find_maximal_match(buff0+seq[i]+j,size-j,t+seq[i],t+seq[i]);
		  t+=h0.len;
		  complexity++;
          if(h0.len>porog) com2num++;
          if(h0.len>porog) com2len+=h0.len;
          }
		  work->add_vector(buff0+seq[i]+j,seq[i]+j,size-j);
		  j++;
		 }
	     while(t<size);
*/
//-----------------	
	work->clear();
	work->set_universe(buff0+seq[i],size);
	for(j=0;j<size;j++)work->add_vector(buff0+seq[i]+j,j,size-j);

	for(j=0;j<howseq;j++)
	   {
	   if(detail_report)fprintf(std_msg,"Decompose %i sequence by %i\n",j,i);
	   if(i!=j) //changed April 2003
	     {
	     t=0;
         com2num=0;
         com2len=0;
		 size1=seq[j+1]-seq[j];
	     complexity=0;
	     while(t<size1)
		  {
		  h0=find_maximal_match(buff0+seq[j]+t,size1-t,size,size);
		  if(detail_report)
		    {
		    h0=determine_sequence(h0);
		    show_item(h0,buff0+seq[j]+t,t);
		    }
		  t+=h0.len;
		  complexity++;
         if(h0.len>porog) com2num++;
         if(h0.len>porog) com2len+=h0.len;
		 }
	     TableDist[j][i]=complexity;//j-th decomposed by i-th
	     TableDist2[j][i]=com2num;
	     TableDist4[j][i]=(double)com2len/(double)size1;
	   }
	   //else complexity=0;//changed - was 1

	   if(detail_report)fprintf(std_msg,"Complexity is %i\n",TableDist[j][i]);
	   fprintf(rep,"%i\t",TableDist[j][i]);
	   }
	fprintf(rep,"\n");
	}
     }
   delete seq;
   }
//Table operations - April 2003

   fprintf(rep,"- Names ------------------ Curr_seq_num=%i\n",Curr_seq_num);
 for(i=0;i<Curr_seq_num;i++)
	{
	for(j=0;j<15;j++)
		 if(Seq_names[i][j]!='\n') fprintf(rep,"%c",Seq_names[i][j]);
	fprintf(rep,"%[%i]\n",seqlengths[i]);
	}
   fprintf(rep,"------------------- \n");

//---------------------
  fprintf(rep,"- TableDist[i][i] - Check --------------------\n");
	for(j=0;j<howseq;j++)
		fprintf(rep,"%i %i\n",j,TableDist[j][j]);

  fprintf(rep,"- Table - Check --------------------\n");
 for(i=0;i<howseq;i++)
	{
	fprintf(rep,"Seq%i(%i)\t",i+1,seqlengths[i]);
	for(j=0;j<howseq;j++)
		fprintf(rep,"%i ",TableDist[i][j]);
	fprintf(rep,"\n");
	}
fprintf(rep,"- Table - Check big component -------------------\n");
 for(i=0;i<howseq;i++)
	{
	fprintf(rep,"Seq%i\t",i+1);
	for(j=0;j<howseq;j++)
		fprintf(rep,"%5i ",TableDist2[i][j]);
	fprintf(rep,"\n");
	}
 fprintf(rep,"- Table - Check proportion of long component -----------------\n");
 for(i=0;i<howseq;i++)
	{
	for(j=0;j<9;j++)
	 if(Seq_names[i][j]!='\n') fprintf(rep,"%c",Seq_names[i][j]);
	 fprintf(rep," ");
	for(j=0;j<howseq;j++)
		fprintf(rep,"%.4f ",TableDist4[i][j]);
	fprintf(rep,"\n");
	}

//---------------------
 for(i=0;i<howseq;i++)
  for(j=0;j<howseq;j++)
	{
    if(TableDist[i][j]>0) TableDistNormalised[i][j]=log((double)seqlengths[i])/(double)TableDist[i][j];
		else TableDistNormalised[i][j]=0;
	}

 for(i=0;i<howseq;i++)
  for(j=0;j<howseq;j++) //was j=i+1
	{
	TableDist2[i][j]=(TableDist[i][j]+TableDist[j][i])/2;
	TableDistNormalised2[i][j]=(TableDistNormalised[i][j]+TableDistNormalised[j][i])/2;
	}

 for(i=0;i<howseq;i++)
	{
   for(j=0;j<howseq;j++) 
	{
	TableDist3[i][j]=TableDist[i][i]+TableDist[j][j]-TableDist[i][j]-TableDist[j][i];
	TableDist4[i][j]=(double)(TableDist3[i][j])/(double)(TableDist[i][i]+TableDist[j][j]);
	}
//	TableDist3[i][i]=0;
//	TableDist4[i][i]=0;
//  a a|b  l(a)/a l(a)/b l(b)/b
//  b|a b
	}
/*
  for(i=0;i<howseq;i++)
  for(j=i;j<howseq;j++) //was j=i+1
	{
	TableDist2[j][i]=TableDist2[i][j];
	TableDistNormalised2[j][i]=TableDistNormalised2[i][j];
	}
*/
//			Seq_names[Curr_seq_num][i]='\n';
 for(i=0;i<Curr_seq_num;i++)
	{
	for(j=0;j<15;j++)
		 if(Seq_names[i][j]!='\n') fprintf(rep,"%c",Seq_names[i][j]);
	fprintf(rep,"%[%i]\n",seqlengths[i]);
	}

  fprintf(rep,"- TableDist --------------------\n");
 for(i=0;i<howseq;i++)
	{
	fprintf(rep,"Seq%i(%i)\t",i+1,seqlengths[i]);
	for(j=0;j<howseq;j++)
		fprintf(rep,"%5i ",TableDist[i][j]);
	fprintf(rep,"\n");
	}
fprintf(rep,"- TableDist2 --------------------\n");
 for(i=0;i<howseq;i++)
	{
	fprintf(rep,"Seq%i\t",i+1);
	for(j=0;j<howseq;j++)
		fprintf(rep,"%5i ",TableDist2[i][j]);
	fprintf(rep,"\n");
	}
fprintf(rep,"- TableDistNormalised --------------------\n");
 for(i=0;i<howseq;i++)
	{
	fprintf(rep,"Seq%i\t",i+1);
	for(j=0;j<howseq;j++)
		fprintf(rep,"%.4f ",TableDistNormalised[i][j]);
	fprintf(rep,"\n");
	}
fprintf(rep,"- TableDistNormalised2 --------------------\n");
 for(i=0;i<howseq;i++)
	{
//	fprintf(rep,"Seq%i\t",i+1);
	for(j=0;j<9;j++)
	 if(Seq_names[i][j]!='\n') fprintf(rep,"%c",Seq_names[i][j]);
	 fprintf(rep," ");
	for(j=0;j<howseq;j++)
		fprintf(rep,"%.4f ",TableDistNormalised2[i][j]);
	fprintf(rep,"\n");
	}

fprintf(rep,"- TableDist3 ([j][j]-[i][j]) --------------------\n");
 for(i=0;i<howseq;i++)
	{
	for(j=0;j<9;j++)
	 if(Seq_names[i][j]!='\n') fprintf(rep,"%c",Seq_names[i][j]);
	 fprintf(rep,"\t");
	for(j=0;j<howseq;j++)
		fprintf(rep,"%i ",TableDist3[i][j]);
	fprintf(rep,"\n");
	}

 fprintf(rep,"- TableDist4 (([j][j]-[i][j])/([i][i]+[j][j]))--------------------\n");
 for(i=0;i<howseq;i++)
	{
	for(j=0;j<9;j++)
	 if(Seq_names[i][j]!='\n') fprintf(rep,"%c",Seq_names[i][j]);
	 fprintf(rep," ");
	for(j=0;j<howseq;j++)
		fprintf(rep,"%.4f ",TableDist4[i][j]);
	fprintf(rep,"\n");
	}
 fprintf(rep,"- TableDist5 (symmetric)--------------------\n");
 for(i=0;i<howseq;i++)
	{
	for(j=0;j<9;j++)
	 if(Seq_names[i][j]!='\n') fprintf(rep,"%c",Seq_names[i][j]);
	 else  fprintf(rep," ");
	 fprintf(rep," ");
	for(j=0;j<howseq;j++)
		fprintf(rep,"%.4f ",(TableDist4[i][j]+TableDist4[j][i])/2);
	fprintf(rep,"\n");
	}

 return 0;
}

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
/*
int open_all_files(FILE **fin,FILE **out,FILE **rep)
{FILE *f=NULL;
 char tempbuff[1024];
 *fin=*out=*rep=NULL;
 if(stdin_flag)
   {
   f=stdin;
   strcpy(main_name,"stdin.txt");
   }else if(main_name[0])f=fopen(main_name,"rb");
 if(f!=NULL)
   {
   *fin=f;
   make_original_name(main_name,res_name,tempbuff);
   strcpy(res_name,tempbuff);
   f=fopen(res_name,"wt");
   if(f!=NULL)
     {
     *out=f;
     if(detail_report)
       {
       make_original_name(main_name,rep_name,tempbuff);
       strcpy(rep_name,tempbuff);
       f=fopen(rep_name,"wt");
       if(f==NULL)return 0;
       *rep=f;
       }else *rep=stdout;
     make_original_name(main_name,second_name,tempbuff);
     strcpy(second_name,tempbuff);
     return 1;
     }else return 0;
   }else return 0;
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


int main(int argc,char *argv[]) //------- console variant !!!
{int i,j;
int flagname=0;
 FILE *fin=NULL,*fsec=NULL,*fout=stdout;

 if(argc>=2)
   {
 
//   use_direct_search=use_inverse_search=1;
//   use_symmetry_search=use_complementary_search=1;
//--- add April - preparation of special statistics --------
for(i=0;i<MAX_N_SPEC_STAT;i++) 
 for(j=0;j<5;j++) 
  special_stat[i][j]=0;
for(i=0;i<MAX_N_SPEC_STAT2;i++) 
 for(j=0;j<6;j++) 
    special_stat_long[i][j]=0;
//July 2003
for(i=0;i<=MAXDISTARRAYSIZE;i++) 
 for(j=0;j<5;j++) 
     {
	 special_dist[i][j]=0;
	 special_dist20[i][j]=0;
     special_dist_sum1[j]=0;
     special_dist_qsum1[j]=0;
     special_dist20_sum1[j]=0;
 	 special_dist20_qsum1[j]=0;
	}

	// analyse command line parameters
    for(i=1;i<argc;i++) //if(handle_key(argv[i])==1)strcpy(main_name,argv[i]);
	{
      switch(handle_key(argv[i]))	//
	    {
	    case 0:	//configuration file
		   break;
				//command line
	    case 1:
			//fprintf(stdout,"Filename:%s (required for console variant)\n",argv[i]);
// This variant for console - April 2002
			if(strstr(main_name,"main_name.seq")!=NULL) 
			 strcpy(main_name,argv[i]);
			else flagname=1; //first name occupied yet

		   if((flagname)&&(strstr(second_name,"second_name.se2")!=NULL) )
			strcpy(second_name,argv[i]);
//		   if((strstr(argv[i],".gif")!=NULL)) gif_filename=argv[i];
//		   seq_filename=argv[i];
		   break;
	    default:fprintf(stdout,"Command line key not recognised. Please, check program parameters!\n");
	    break;
		}
	 }
   use_forward=(use_direct_search||use_complementary_search);
   use_reverse=(use_inverse_search||use_symmetry_search);

if(!(use_forward||use_reverse) ) 
{
fprintf(std_msg,"Repeat search not defined! Use at least one type of repeats: direct,inverted, etc.\n");
return 0;
}

//   if(open_all_files(&fin,&fout,&std_msg))
 fin=fopen(main_name,"rb");
if(count_method==2) 
 fsec=fopen(second_name,"rb");
// else if(!show_help)fprintf(std_msg,"Filename required or can't open files\n");

     {
     buff0=new byte[SEQUENCE_LEN+2];
     if(buff0!=NULL)
       {
       base=new genfileinput();
       base->input_init(alphabet_pointer);
       if(base->success_init())
	 {

		ALPHABET_LEN=base->get_alphabet_len();
	 if((!best_function)&&(count_method!=4)&&(deffunc[0]))
	   {
	   if(!set_complementary_function())fprintf(std_msg,"Bad complementary function\n");
	   }
	 if(best_function) work=new multidivizion(NSIZE,ALPHABET_LEN);
	 else
		{
		if((count_method==0)&&(SLIDE_STEP==1))work=new slide_divizion(NSIZE,ALPHABET_LEN);
		  else work=new divizion(NSIZE,ALPHABET_LEN);
		}

	 if(work->success_init())
	   {
	     if(seq_init())
	       {
	     switch(count_method) //main calculation of complexity profile
		   {
		   case 0:calculate_profile(fin,fout);
			  break;
		   case 1:calculate_whole(fin,fout);
			  break;
		   case 2:buff1=new byte[SEQUENCE_LEN+2];
			  if(buff1!=NULL)
			    {
			    fsec=fopen(second_name,"rb");
			    if(fsec!=NULL)
			      {
			      calculate_analyse(fsec,fin,fout);
			      fclose(fsec);
			      }else fprintf(std_msg,"File open error\n");
			    delete buff1;
			    }else fprintf(std_msg,"Not enought memory for second buffer\n");
			  break;
		   case 3:
//			    calc_whole_seqs(fin,fout);
//                fseek(fin,0,SEEK_SET);
				calculate_table(fin,fout);
			  break;
		   case 4:if(best_function)fprintf(std_msg,"Can't count vector complexity with best function\n");
			    else calculate_vector(fin,fout);
			  break;
//		   case 5:calculateFASTAprofile(fin,fout);
//			  break;
		   default:
			   fprintf(std_msg,"Method is not defined/not recognised from command line\n");
            break;
		 }//switch
	      }else fprintf(std_msg,"Not enought memory for decompose subsystem\n");

		}else fprintf(std_msg,"Not enought memory for decompose subsystem\n");

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





//-------------
