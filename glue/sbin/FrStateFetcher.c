/***************************
 Example Usage:
 FrStateFeetcher file1 file2 ... filen

 prints (for each file)
 <startgps> <startgps ns> <end gps> <end gps ns> <statevec val> <statevec mantissa> <segnum> <ifo> <filename>

 where the <..gps> + <...gps ns> == true same start/end gps time.
 With the usual double-precision floating point precision problems.

 Originally written by,
 Ben Johnson
 LIGO Hanford Observatory
 johnson_b@spammers.will.never.find.me.ligo-wa.caltech.edu


********************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <FrameL.h>

/* default string size */
#define STRING_SIZE 1024

/* Name of channels */
#define H1_SV "H1:IFO-SV_STATE_VECTOR"
#define H2_SV "H2:IFO-SV_STATE_VECTOR"
#define L1_SV "L1:IFO-SV_STATE_VECTOR"
#define G1_SV "G1:DER_DATA_QUALITY"
#define V1_SV "V1:Hrec_veto_dataQuality"

#define H1_SEG "H1:IFO-SV_SEGNUM"
#define H2_SEG "H2:IFO-SV_SEGNUM"
#define L1_SEG "L1:IFO-SV_SEGNUM"


/* max number of interferometers that can be
 processed by this code */
#define IFO_NUMBER  512

/* Maximum number of recursions to allow segment_number
   mode code to run. */
#define MAX_RECUSIONS 1024

typedef char* IfoList;
/* holds state information for
 *    segment recursion structure */
typedef struct tagsegnum_struct {
    FRULONG nelem;
    float* vals;
    FRULONG* count;
} segnum_struct;



/* For segment number mode calculations */
float segnum_mode(FrVect*, unsigned long, unsigned long);
float grab_mode(segnum_struct*);
int init_segnum_struct(FRULONG,segnum_struct*);
void del_segnum_struct(segnum_struct*);
void increment_mode(float,segnum_struct*);


/* ifoMap header, see function body for details */
int ifoMap(IfoList*,char*);

int Print_State_H1(char* filename,FILE* fp);
int Print_State_H2(char* filename,FILE* fp);
int Print_State_L1(char* filename,FILE* fp);
int Print_State_G1(char* filename,FILE* fp);

/* Usage */
void print_usage(void);

/* Because my Solaris box doesn't have round() */
double myround(double);

/* BEGIN MAIN */
int main(int argc,char**argv)
{
    size_t i,k;
    unsigned int science,common;
    char filename[STRING_SIZE];
    char *detname;
    FrFile *iFile;
    FrTOC* toc;
    FrTOCdet *det;
    IfoList* ifo;
  

    /* command line argument check */
    if(argc < 2)
    {
        fprintf(stderr,"You must supply at least one file name.\n\n");
        print_usage();
        exit(1);
    }

    if( strncmp("-h",argv[1],2) == 0 )
    {
        print_usage();
        exit(0);
    }


    /* loop over each file */
    for(i=1;i<argc;i++)
    {

        strncpy(filename,argv[i],STRING_SIZE-1);
     
        /* open file, make sure it exists */
        iFile = FrFileINew(filename);
    
        if(iFile == NULL)
        {
            fprintf(stderr,"%s: %s: Could not open file\n",argv[0],filename);
            continue;
        }

        ifo = (IfoList*)calloc(IFO_NUMBER,sizeof(IfoList));
        toc = FrTOCReadFull(iFile);

        if(toc == NULL)
        {
            fprintf(stderr,"%s: %s: File does not contain toc\n",argv[0],filename);
            FrFileIClose(iFile); /* frees the toc as well ? */
            free(ifo);
            continue;
        }
        if(toc->detector == NULL)
        {
            fprintf(stderr,"%s: %s: File does not contain \"detector\" field in toc\n",argv[0],filename);
            FrFileIClose(iFile);   /* frees the toc as well ? */
            free(ifo);
            continue;
        }

        det = toc->detector;
        while(det != NULL)
        {
            detname = det->name;
            ifoMap(ifo,detname);
            det = det->next;
        }
        FrFileIClose(iFile); /* frees the toc as well ? */
        science = 0;
        common = 0;
        for(k=0; ifo[k] != NULL; k++)
        {

            if( strncmp(ifo[k],"H1",STRING_SIZE) == 0 )
            {
                /* Get H1 data */
                Print_State_H1(filename,stdout);
            }
            else if( strncmp(ifo[k],"H2",STRING_SIZE) == 0 )
            {
                Print_State_H2(filename,stdout);
            }
            else if( strncmp(ifo[k],"L1",STRING_SIZE) == 0 )
            {
                Print_State_L1(filename,stdout);
            }
            else if( strncmp(ifo[k],"G1",STRING_SIZE) == 0 )
            {
                Print_State_G1(filename,stdout);
            }
	    else if( strncmp(ifo[k], "V1", STRING_SIZE) == 0 )
	    {
		 Print_State_V1(filename,stdout);
	    }
            
            free(ifo[k]);
        }
        free(ifo);

    }

    exit(0);
}

void print_usage(void)
{
    fprintf(stdout,"The program prints state vector and segment number channels in given file list.\n");
    fprintf(stdout,"The channels examined are, Xn:IFO-SV_STATE_VECTOR, and\n");
    fprintf(stdout,"Xn:IFO-SV_SEGNUM, for the state vector and segment number\n");
    fprintf(stdout,"respectively. (Xn, can be H1, H2, or L1 currently)\n");
    fprintf(stdout,"\nThe columns are as follows:\n");
    fprintf(stdout,"Col 1) Gps start of state vector segment.\n");
    fprintf(stdout,"Col 2) Gps start nanoseconds.\n");
    fprintf(stdout,"Col 3) Gps end of state vector segment.\n");
    fprintf(stdout,"Col 4) Gps end nanoseconds.\n");
    fprintf(stdout,"Col 5) State vector value (integer part).\n");
    fprintf(stdout,"Col 6) State vector value (fractional part).\n");
    fprintf(stdout,"Col 7) Mode of the segment number during printed GPS interval.\n");
    fprintf(stdout,"Col 8) Interferometer.\n");
    fprintf(stdout,"Col 9) Filename from where information was pulled.\n\n");

}


int ifoMap(IfoList*ifo,char*name)
{
    /* Given name, maps the name to the appropriate
     interferometer name, and if the entry does not
     already exist, adds it to the IfoList.

     mapping is stored in MapDict* dict;

     This function returns the number of entries added
     to the IfoList. ( i.e. 0 or 1)
     */

    size_t i;
    int modh1,modh2;

    /* GEO, G1 */
    if( ! strncmp("GEO",name,STRING_SIZE) )
    {
        for(i=0; ifo[i] != NULL; i++)
        {
            if( ! strncmp("G1",ifo[i],STRING_SIZE))
            {
                return 0;
            }
        }
        ifo[i] = (char*)malloc(STRING_SIZE*sizeof(char));
        strncpy(ifo[i],"G1",STRING_SIZE);
        return 1;
    }

    /* LHO, H1, H2 --> Due to old (i.e. many before 2004) frames */
    if( ! strncmp("LHO",name,STRING_SIZE) )
    {
        modh1 = modh2 = 1;
        for(i=0; ifo[i] != NULL; i++)
        {
            if( ! strncmp("H1",ifo[i],STRING_SIZE) )
            {
                modh1 = 0;
            }
            if( ! strncmp("H2",ifo[i],STRING_SIZE) )
            {
                modh2 = 0;
            }
        }
        if(modh1)
        {
            ifo[i] = (char*)malloc(STRING_SIZE*sizeof(char));
            strncpy(ifo[i],"H1",STRING_SIZE);
            i++;
        }
        if(modh2)
        {
            ifo[i] = (char*)malloc(STRING_SIZE*sizeof(char));
            strncpy(ifo[i],"H2",STRING_SIZE);
            i++;
        }

        return modh1 + modh2;
    }

    /* LHO_4k, H1 */
    if( ! strncmp("LHO_4k",name,STRING_SIZE))
    {
        for(i=0; ifo[i] != NULL; i++)
        {
            if( ! strncmp("H1",ifo[i],STRING_SIZE))
            {
                return 0;
            }
        }
        ifo[i] = (char*)malloc(STRING_SIZE*sizeof(char));
        strncpy(ifo[i],"H1",STRING_SIZE);
        return 1;
    }

    /* LHO_2k, H2 */
    if( ! strncmp("LHO_2k",name,STRING_SIZE))
    {
        for(i=0; ifo[i] != NULL; i++)
        {
            if( ! strncmp("H2",ifo[i],STRING_SIZE))
            {
                return 0;
            }
        }
        ifo[i] = (char*)malloc(STRING_SIZE*sizeof(char));
        strncpy(ifo[i],"H2",STRING_SIZE);
        return 1;
    }

    /* LLO_4k, L1 */
    if( ! strncmp("LLO_4k",name,STRING_SIZE))
    {
        for(i=0; ifo[i] != NULL; i++)
        {
            if( ! strncmp("L1",ifo[i],STRING_SIZE))
            {
                return 0;
            }
        }
        ifo[i] = (char*)malloc(STRING_SIZE*sizeof(char));
        strncpy(ifo[i],"L1",STRING_SIZE);
        return 1;
    }

    /* Virgo, V1, or V2 ???? */
    if( ! strncmp("Virgo",name,STRING_SIZE))
    {
        for(i=0; ifo[i] != NULL; i++)
        {
            if( ! strncmp("V1",ifo[i],STRING_SIZE))
            {
                return 0;
            }
        }
        ifo[i] = (char*)malloc(STRING_SIZE*sizeof(char));
        strncpy(ifo[i],"V1",STRING_SIZE);
        return 1;
    }

    return 0;

}

int Print_State_G1(char* filename,FILE* fp)
{
  
    FrFile *iFile;
    FrVect *vect;
    FRULONG j,nData;
    unsigned short sample,start_sample=0;
    unsigned short *ptr;
    double start_time=0.0,end_time=0.0;
    double file_init, file_final;
    char *ifo_name = "G1";
    char sStart[22],newS[10];
    char sEnd[22],newE[10];
    size_t len,i,k;
    
    iFile = FrFileINew(filename);
    
    file_init = FrFileITStart(iFile);
    file_final = FrFileITEnd(iFile);

    vect = FrFileIGetVect(iFile,G1_SV,file_init,file_final - file_init);

    while(vect != NULL)
    {
        nData = vect->nData;

        ptr = (unsigned short*)vect->data;
        start_sample = ptr[0];
        start_time = vect->GTime;
        end_time = start_time + vect->dx[0];

        for(j=1; j < nData; j++)
        {

            sample = ptr[j];

            if(sample != start_sample)
            {

                /* formatting for nanoseconds field */
                snprintf(sStart,21,"%0.9lf",start_time);
                snprintf(sEnd,21,"%0.9lf",end_time);

                len = strlen(sStart);
                for(i = len,k=9; i > len - 10 && k >=0; i--,k--)
                {
                    newS[k] = sStart[i];
                }

                len = strlen(sEnd);
                for(i = len,k=9; i > len - 10 && k >=0; i--,k--)
                {
                    newE[k] = sEnd[i];
                }


                fprintf(fp,"%0.lf %s %0.lf %s %d %d %d %s %s\n",start_time,newS,end_time,
                        newE,start_sample,(int)0,(int)-1,ifo_name,filename);


                start_sample = sample;
                start_time = end_time;
            }
            end_time += vect->dx[0];
        }
        vect = vect->next;
    }


    /* formatting for nanoseconds field */
    snprintf(sStart,21,"%0.9lf",start_time);
    snprintf(sEnd,21,"%0.9lf",end_time);


    len = strlen(sStart);
    for(i = len,k=9; i > len - 10 && k >=0; i--,k--)
    {
        newS[k] = sStart[i];
    }

    len = strlen(sEnd);
    for(i = len,k=9; i > len - 10 && k >=0; i--,k--)
    {
        newE[k] = sEnd[i];
    }

    fprintf(fp,"%0.lf %s %0.lf %s %d %d %d %s %s\n",start_time,newS,end_time,
            newE,start_sample,(int)0,(int)-1,ifo_name,filename);

    FrVectFree(vect);
    
    FrFileIEnd(iFile);
    return 0;
   
}

int Print_State_V1(char* filename,FILE* fp)
{
  
    FrFile *iFile;
    FrVect *vect;
    FRULONG j,nData;
    float sample,start_sample=0;
    float *ptr;
    double start_time=0.0,end_time=0.0;
    double file_init, file_final;
    char *ifo_name = "V1";
    char sStart[22],newS[10];
    char sEnd[22],newE[10];
    size_t len,i,k;
    
    iFile = FrFileINew(filename);
    
    file_init = FrFileITStart(iFile);
    file_final = FrFileITEnd(iFile);

    vect = FrFileIGetVect(iFile,V1_SV,file_init,file_final - file_init);

    while(vect != NULL)
    {
        nData = vect->nData;

        ptr = (float*)vect->data;
        start_sample = ptr[0];
        start_time = vect->GTime;
        end_time = start_time + vect->dx[0];

        for(j=1; j < nData; j++)
        {

            sample = ptr[j];

            if(sample != start_sample)
            {

                /* formatting for nanoseconds field */
                snprintf(sStart,21,"%0.9lf",start_time);
                snprintf(sEnd,21,"%0.9lf",end_time);

                len = strlen(sStart);
                for(i = len,k=9; i > len - 10 && k >=0; i--,k--)
                {
                    newS[k] = sStart[i];
                }

                len = strlen(sEnd);
                for(i = len,k=9; i > len - 10 && k >=0; i--,k--)
                {
                    newE[k] = sEnd[i];
                }


                fprintf(fp,"%0.lf %s %0.lf %s %0.f %d %d %s %s\n",start_time,newS,end_time,
                        newE,start_sample,(int)0,(int)-1,ifo_name,filename);


                start_sample = sample;
                start_time = end_time;
            }
            end_time += vect->dx[0];
        }
        vect = vect->next;
    }


    /* formatting for nanoseconds field */
    snprintf(sStart,21,"%0.9lf",start_time);
    snprintf(sEnd,21,"%0.9lf",end_time);


    len = strlen(sStart);
    for(i = len,k=9; i > len - 10 && k >=0; i--,k--)
    {
        newS[k] = sStart[i];
    }

    len = strlen(sEnd);
    for(i = len,k=9; i > len - 10 && k >=0; i--,k--)
    {
        newE[k] = sEnd[i];
    }

    fprintf(fp,"%0.lf %s %0.lf %s %0.f %d %d %s %s\n",start_time,newS,end_time,
            newE,start_sample,(int)0,(int)-1,ifo_name,filename);

    FrVectFree(vect);
    
    FrFileIEnd(iFile);
    return 0;
   
}

int Print_State_H1(char* filename,FILE*fp)
{
    FrFile *iFile;
    FrAdcData *adc,*seg;
    FrVect *vect,*segvect;
    FRULONG j,nData,sidx;
    float sample,start_sample,sn_mode;
    float *ptr;
    double start_time,end_time;
    char *ifo_name = "H1";
    char sStart[22],newS[10];
    char sEnd[22],newE[10];
    char sSamp[22],newSamp[10];
    size_t len,i,k;

  
    iFile = FrFileINew(filename);
    adc = FrAdcDataReadT(iFile,H1_SV,0);
    seg = FrAdcDataReadT(iFile,H1_SEG,0);

    if(adc != NULL && seg != NULL)
    {
        vect = adc->data;
        segvect = seg->data;

        /* Current algorithm requires the state vector and
         segnum channels to have the same number of samples */
        if(vect-> nData != segvect -> nData)
        {
            fprintf(stderr,"The channels, %s and %s, do not have the same sample rate.\n\
                    This code is unable to produce correct output when the sample\n\
                    rate of those two channels differs. Found in file %s.\n",H1_SV,H1_SEG,filename);
            FrAdcDataFree(adc);
            FrAdcDataFree(seg);
            return 1;
        }

        while(vect != NULL)
        {
            nData = vect->nData;
            
            sidx = 0;
            ptr = (float*)vect->data;
            start_sample = ptr[0];
            start_time = vect->GTime;

            end_time = start_time + vect->dx[0];

            for(j=1; j < nData; j++)
            {

                sample = ptr[j];

                if(sample != start_sample)
                {


                    /* formatting for nanoseconds field */
                    snprintf(sStart,21,"%0.9lf",start_time);
                    snprintf(sEnd,21,"%0.9lf",end_time);

                    /* formatting for state vector mantissa */
                    snprintf(sSamp,21,"%0.9f",start_sample);

                    len = strlen(sStart);
                    for(i = len,k=9; i > len - 10 && k >=0; i--,k--)
                    {
                        newS[k] = sStart[i];
                    }

                    len = strlen(sEnd);
                    for(i = len,k=9; i > len - 10 && k >=0; i--,k--)
                    {
                        newE[k] = sEnd[i];
                    }

                    len = strlen(sSamp);
                    for(i = len,k=9; i > len - 10 && k >=0; i--,k--)
                    {
                        newSamp[k] = sSamp[i];
                    }
                    
                    /* Find mode of segment number */
                    sn_mode = segnum_mode(segvect,sidx,j);
                    fprintf(fp,"%0.lf %s %0.lf %s %0.lf %s %0.f %s %s\n",start_time,newS,end_time,
                            newE,start_sample,newSamp,sn_mode,ifo_name,filename);


                    start_sample = sample;
                    start_time = end_time;
                    sidx = j;
                }
                end_time += vect->dx[0];
            }
            vect = vect->next;
        }


        /* formatting for nanoseconds field */
        snprintf(sStart,21,"%0.9lf",start_time);
        snprintf(sEnd,21,"%0.9lf",end_time);

        /* formatting for state vector mantissa */
        snprintf(sSamp,21,"%0.9f",start_sample);

        len = strlen(sStart);
        for(i = len,k=9; i > len - 10 && k >=0; i--,k--)
        {
            newS[k] = sStart[i];
        }

        len = strlen(sEnd);
        for(i = len,k=9; i > len - 10 && k >=0; i--,k--)
        {
            newE[k] = sEnd[i];
        }

        len = strlen(sSamp);
        for(i = len,k=9; i > len - 10 && k >=0; i--,k--)
        {
            newSamp[k] = sSamp[i];
        }
        /* Find mode of segment number */
        sn_mode = segnum_mode(segvect,sidx,j);

        fprintf(fp,"%0.lf %s %0.lf %s %0.lf %s %0.f %s %s\n",start_time,newS,end_time,
                newE,start_sample,newSamp,sn_mode,ifo_name,filename);

        FrAdcDataFree(adc);
        FrAdcDataFree(seg);
    }
    FrFileIEnd(iFile);
    return 0;
}


int Print_State_H2(char* filename,FILE*fp)
{
    FrFile *iFile;
    FrAdcData *adc,*seg;
    FrVect *vect,*segvect;
    FRULONG j,nData,sidx;
    float sample,start_sample,sn_mode;
    float *ptr;
    double start_time,end_time;
    char *ifo_name = "H2";
    char sStart[22],newS[10];
    char sEnd[22],newE[10];
    char sSamp[22],newSamp[10];
    size_t len,i,k;

  
    iFile = FrFileINew(filename);
    adc = FrAdcDataReadT(iFile,H2_SV,0);
    seg = FrAdcDataReadT(iFile,H2_SEG,0);

    if(adc != NULL && seg != NULL)
    {
        vect = adc->data;
        segvect = seg->data;

        /* Current algorithm requires the state vector and
        segnum channels to have the same number of samples */
        if(vect-> nData != segvect -> nData)
        {
            fprintf(stderr,"The channels, %s and %s, do not have the same sample rate.\n\
                    This code is unable to produce correct output when the sample\n\
                            rate of those two channels differs. Found in file %s.\n",H1_SV,H1_SEG,filename);
            FrAdcDataFree(adc);
            FrAdcDataFree(seg);
            return 1;
        }

        while(vect != NULL)
        {
            nData = vect->nData;
            
            sidx = 0;
            ptr = (float*)vect->data;
            start_sample = ptr[0];
            start_time = vect->GTime;

            end_time = start_time + vect->dx[0];

            for(j=1; j < nData; j++)
            {

                sample = ptr[j];

                if(sample != start_sample)
                {


                    /* formatting for nanoseconds field */
                    snprintf(sStart,21,"%0.9lf",start_time);
                    snprintf(sEnd,21,"%0.9lf",end_time);

                    /* formatting for state vector mantissa */
                    snprintf(sSamp,21,"%0.9f",start_sample);

                    len = strlen(sStart);
                    for(i = len,k=9; i > len - 10 && k >=0; i--,k--)
                    {
                        newS[k] = sStart[i];
                    }

                    len = strlen(sEnd);
                    for(i = len,k=9; i > len - 10 && k >=0; i--,k--)
                    {
                        newE[k] = sEnd[i];
                    }

                    len = strlen(sSamp);
                    for(i = len,k=9; i > len - 10 && k >=0; i--,k--)
                    {
                        newSamp[k] = sSamp[i];
                    }
                    
                    /* Find mode of segment number */
                    sn_mode = segnum_mode(segvect,sidx,j);
                    fprintf(fp,"%0.lf %s %0.lf %s %0.lf %s %0.f %s %s\n",start_time,newS,end_time,
                            newE,start_sample,newSamp,sn_mode,ifo_name,filename);


                    start_sample = sample;
                    start_time = end_time;
                    sidx = j;
                }
                end_time += vect->dx[0];
            }
            vect = vect->next;
        }


        /* formatting for nanoseconds field */
        snprintf(sStart,21,"%0.9lf",start_time);
        snprintf(sEnd,21,"%0.9lf",end_time);

        /* formatting for state vector mantissa */
        snprintf(sSamp,21,"%0.9f",start_sample);

        len = strlen(sStart);
        for(i = len,k=9; i > len - 10 && k >=0; i--,k--)
        {
            newS[k] = sStart[i];
        }

        len = strlen(sEnd);
        for(i = len,k=9; i > len - 10 && k >=0; i--,k--)
        {
            newE[k] = sEnd[i];
        }

        len = strlen(sSamp);
        for(i = len,k=9; i > len - 10 && k >=0; i--,k--)
        {
            newSamp[k] = sSamp[i];
        }
        /* Find mode of segment number */
        sn_mode = segnum_mode(segvect,sidx,j);

        fprintf(fp,"%0.lf %s %0.lf %s %0.lf %s %0.f %s %s\n",start_time,newS,end_time,
                newE,start_sample,newSamp,sn_mode,ifo_name,filename);

        FrAdcDataFree(adc);
        FrAdcDataFree(seg);
    }
    FrFileIEnd(iFile);
    return 0;
}

int Print_State_L1(char* filename,FILE*fp)
{
    FrFile *iFile;
    FrAdcData *adc,*seg;
    FrVect *vect,*segvect;
    FRULONG j,nData,sidx;
    float sample,start_sample,sn_mode;
    float *ptr;
    double start_time,end_time;
    char *ifo_name = "L1";
    char sStart[22],newS[10];
    char sEnd[22],newE[10];
    char sSamp[22],newSamp[10];
    size_t len,i,k;

  
    iFile = FrFileINew(filename);
    adc = FrAdcDataReadT(iFile,L1_SV,0);
    seg = FrAdcDataReadT(iFile,L1_SEG,0);

    if(adc != NULL && seg != NULL)
    {
        vect = adc->data;
        segvect = seg->data;

        /* Current algorithm requires the state vector and
        segnum channels to have the same number of samples */
        if(vect-> nData != segvect -> nData)
        {
            fprintf(stderr,"The channels, %s and %s, do not have the same sample rate.\n\
                    This code is unable to produce correct output when the sample\n\
                            rate of those two channels differs. Found in file %s.\n",H1_SV,H1_SEG,filename);
            FrAdcDataFree(adc);
            FrAdcDataFree(seg);
            return 1;
        }

        while(vect != NULL)
        {
            nData = vect->nData;
            
            sidx = 0;
            ptr = (float*)vect->data;
            start_sample = ptr[0];
            start_time = vect->GTime;

            end_time = start_time + vect->dx[0];

            for(j=1; j < nData; j++)
            {

                sample = ptr[j];

                if(sample != start_sample)
                {


                    /* formatting for nanoseconds field */
                    snprintf(sStart,21,"%0.9lf",start_time);
                    snprintf(sEnd,21,"%0.9lf",end_time);

                    /* formatting for state vector mantissa */
                    snprintf(sSamp,21,"%0.9f",start_sample);

                    len = strlen(sStart);
                    for(i = len,k=9; i > len - 10 && k >=0; i--,k--)
                    {
                        newS[k] = sStart[i];
                    }

                    len = strlen(sEnd);
                    for(i = len,k=9; i > len - 10 && k >=0; i--,k--)
                    {
                        newE[k] = sEnd[i];
                    }

                    len = strlen(sSamp);
                    for(i = len,k=9; i > len - 10 && k >=0; i--,k--)
                    {
                        newSamp[k] = sSamp[i];
                    }
                    
                    /* Find mode of segment number */
                    sn_mode = segnum_mode(segvect,sidx,j);
                    fprintf(fp,"%0.lf %s %0.lf %s %0.lf %s %0.f %s %s\n",start_time,newS,end_time,
                            newE,start_sample,newSamp,sn_mode,ifo_name,filename);


                    start_sample = sample;
                    start_time = end_time;
                    sidx = j;
                }
                end_time += vect->dx[0];
            }
            vect = vect->next;
        }


        /* formatting for nanoseconds field */
        snprintf(sStart,21,"%0.9lf",start_time);
        snprintf(sEnd,21,"%0.9lf",end_time);

        /* formatting for state vector mantissa */
        snprintf(sSamp,21,"%0.9f",start_sample);

        len = strlen(sStart);
        for(i = len,k=9; i > len - 10 && k >=0; i--,k--)
        {
            newS[k] = sStart[i];
        }

        len = strlen(sEnd);
        for(i = len,k=9; i > len - 10 && k >=0; i--,k--)
        {
            newE[k] = sEnd[i];
        }

        len = strlen(sSamp);
        for(i = len,k=9; i > len - 10 && k >=0; i--,k--)
        {
            newSamp[k] = sSamp[i];
        }
        /* Find mode of segment number */
        sn_mode = segnum_mode(segvect,sidx,j);

        fprintf(fp,"%0.lf %s %0.lf %s %0.lf %s %0.f %s %s\n",start_time,newS,end_time,
                newE,start_sample,newSamp,sn_mode,ifo_name,filename);

        FrAdcDataFree(adc);
        FrAdcDataFree(seg);
    }
    FrFileIEnd(iFile);
    return 0;
}

float segnum_mode(FrVect* segnum_channel, FRULONG startidx, FRULONG stopidx)
{
    /* This function is designed to find and return
    the mode of the segnum channel between the closed 
    interval defined by startidx and stopidx.
        
    The algorithm is designed to take advantage of the
    segment number's propensity to take on, at most, two
    different values in the given interval. The algorithm
    will work in the general case, but not be very efficient
    at it.
    
    */
   
    float answer;
    int i;
    segnum_struct counters;
    float* dF;
    
    /* die if start and end indices are not sane */
    if(startidx > stopidx) {
        return (float)-2.0; 
    }
    if(stopidx > segnum_channel->nData) {
        return (float)-3.0;
    }
    
    if( init_segnum_struct(segnum_channel->nData, &counters) < 0)
    {
        perror("init_segnum_struct() failed");
        exit(1);
    }
    
    dF = (float*)segnum_channel->data;
    for(i = startidx; i < stopidx; i++)
    {
        increment_mode(dF[i],&counters);
    }
    
    answer = grab_mode(&counters);
    del_segnum_struct(&counters);
    
    return answer;
}

int init_segnum_struct(FRULONG nData,segnum_struct* counters)
{
    counters->nelem = 0;
    counters->vals = (float*)malloc(nData * sizeof(float));
   
    if(counters->vals == NULL)
    {
        return -1; 
    }
    counters->count = (FRULONG*)malloc(nData * sizeof(FRULONG));
    if(counters->count == NULL)
    {
        return -2;
    }
    return 0;
}

void del_segnum_struct(segnum_struct* n)
{
    free(n->vals);
    free(n->count);
}

void increment_mode(float newval,segnum_struct* n)
{
    FRULONG i;
    
    for(i=0;i<n->nelem;i++)
    {
        if(n->vals[i] == newval)
        {
            (n->count[i])++;
            return;
        }
    }
   
    n->vals[i] = newval;
    n->count[i] = 1;
    (n->nelem)++;
}

float grab_mode(segnum_struct* n)
{
    FRULONG i;
    FRULONG cur_max = 0;
    float answer = -1234;
    
    for(i=0; i < n->nelem; i++)
    {
        
        if(n->count[i] >= cur_max)
        {
            cur_max = n->count[i];
            answer = n->vals[i];
           
        }
    }
    
    
    return answer;
}


/*
 I'm sure this has the 5 is not a power of 2
 problem... but since my Solaris box doesn't
 have round()
 */
double myround(double x)
{
    double flr_x;

    flr_x = floor(x);

    if(x + 0.5 >= flr_x + 1.0)
    {
        return flr_x + 1.0;
    } else {
        return flr_x;
    }
}
