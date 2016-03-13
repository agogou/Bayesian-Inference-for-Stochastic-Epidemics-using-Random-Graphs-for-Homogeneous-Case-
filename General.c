//
//  main.c
//  SIR_graph_inference
//
//

#include<stdio.h>    /* printf */
#include<stdlib.h>   /* abs */
#include<math.h>
#include<time.h>
#include <stdbool.h>

int look2(int*vector,int**matrix,int n,int j);
void award(int*A,int**B,int*A2,int**B2,int n);
int check(int **c,int *d,int start,int n);
void dohits(int **c,int *d,int start,int n,int *ht);
long int factorial(int a);
int sum(int *pa, int num_elements);
long int choose(int n,int k);
double randn();
double randomInteger01(void);
int randomIntegern(int n);
void loglikelihood(double *result, double beta,int **c,int *d,int n,int pop);
void display_mat_deg(int ** c, int * d, int n);

int main()
{
    int n=30,pop=120;
    int i,j,k,l,pap,*d,*propd;
    double sbeta=0.01,beta=0.001,mean_beta=0,R0=0;//50=0.004 25=0.007
    double prop,numer[2],denom[2],u,acc,sumll;
    long double pap2;
    int burn=10000,samgap=10000,samsize=50,looplen,counter,p,miko=0,run=0;
    int A=0,B=0,A2=0,B2=0;
    bool done;
    int **c, **propc;
    FILE *fd;
    
    /*dynamic allocation of memory for c*/
    c= (int **)malloc(n * sizeof(int *)); 
    for (i=0; i<n; i++)
    {
        c[i] = (int *)malloc(n * sizeof(int));
    }
    if(c== NULL)
    {
        printf("Not enough memory\n");
        exit(1);
    }
   
    propc=(int **)malloc(n * sizeof(int *));  
    for (i=0; i<n; i++)
    {
        propc[i] = (int *)malloc(n * sizeof(int));
    }
    if(propc== NULL)
    {
        printf("Not enough memory\n");
        exit(1);
    }
    d=(int *) malloc(n*sizeof(int));
    propd=(int *) malloc(n*sizeof(int));
    
    /*fill the c and d matrices with 0's*/
    for(i=0; i<n; i++)      
    {
        for(j=0; j<n; j++)
            c[i][j]=0;
        d[i]=0;
   	}
    
    /*define the initial c matrix by a chain from individual 1 until the n, all node have degree 1*/
    for(i=0; i<n-1; i++)   
    {
        c[i][0]=i+2;
        d[i]=1;
    }
    
    display_mat_deg(c, d, n);
    
    /*initialization of the random generator using the internal clock*/
    srand((unsigned int)time(NULL));
    srand(57);
    
    /*opens the file for saving results*/
    if((fd=fopen("arxeio30.txt","w"))==NULL)  
        printf("File could not opened\n"); 
    
    burn=burn-(2*burn);
    looplen=samgap*samsize;
    
    for(counter=burn; counter<=looplen; counter++) //start
    {        
        prop=randn(beta,sbeta); //draws a proposal for beta
        prop=fabs(prop);        //takes the absolute value
        loglikelihood(numer,prop,c,d,n,pop);
        if(numer[1]==1)
        {
            loglikelihood(denom,beta,c,d,n,pop);
            acc=exp(numer[0]-denom[0]);
            u=randomInteger01();
            if(u<acc)
            {
                miko++;
                beta=prop;
            }
        }
        
        p=randomIntegern(2);                  //if 1->add link if 2->delete link  
        
        if((p==1) && (sum(d,n)<n*(n-1)))    
        {   A++;
            do
            {
                done=true;
                k=randomIntegern(n);   
                l=randomIntegern(n);   
                if((l==k) || (d[k-1]==n-1))  
                {
                    done=false;    
                }
                else if(d[k-1]>=1)       /* we control if already exist */
                {                                          
                    for(i=0; i<d[k-1]; i++)
                    {
                        if(c[k-1][i]==l)
                        {
                            done=false;
                        }
                    }
                }
            }while(done==false);

            award(propd,propc,d,c,n);
            //printf("\nwe choose the %d individual to ADD the %d link\n",k,l);
            //look2(propd,propc,n,k);
            propd[k-1]++;         
            propc[k-1][d[k-1]]=l;  /* we dreaw the new link */
            //printf("egine-->\n");
            //look2(propd,propc,n,k);            
       
            loglikelihood(numer,beta,propc,propd,n,pop);
            //printf("ADD link valid: %f\n",numer[1]);
            if(numer[1]==1)
            {  
                acc=(exp(beta)-1)*(n*(n-1)-sum(d,n))/(sum(d,n)+1);
                u=randomInteger01();
                //printf("Add u=%f ,acc= %f \n",u,acc);
                if(u<acc)
                {   A2++;
                    award(d,c,propd,propc,n);  
                }
            }
        }
        else if(p==2 && sum(d,n)>(n-1)) /* if 2 delete a link */
        {   B++;
            i=randomIntegern(sum(d,n));
            k=1;
            l=0;
            while(l+d[k-1]<i)
            {
                l=l+d[k-1];
                k++;
            }
            l=i-l;
            
            award(propd,propc,d,c,n);
            //printf("\nwe choose for %d individual to DELETE the %d link",k,l);
            //look2(propd,propc,n,k); 
            //printf("egine-->\n");
            if(l<d[k-1] && l>0)
            {
                for(i=(l-1); i<=(d[k-1]-1); i++)
                {
                    propc[k-1][i]=propc[k-1][i+1]; 
                }
            }
            propc[k-1][d[k-1]-1]=0;
            propd[k-1]--;
            
            //look2(propd,propc,n,k); 
            //printf("(%d,%d) -X->\n",k,l);
            //display_mat_deg(propc,propd, n);
            
            loglikelihood(numer,beta,propc,propd,n,pop);
            //printf("DELETE LINK valid: %f\n",numer[1]);
            if(numer[1]==1)  
            {   
                loglikelihood(denom,beta,c,d,n,pop);
                acc=exp(numer[0]-denom[0]);
                acc=acc*sum(d,n)/(n*(n-1)-sum(d,n)+1);
                u=randomInteger01();
                //printf("delete: u=%f ,acc= %f \n",u,acc);
                if(u<acc)
                {    B2++;
                     award(d,c,propd,propc,n); 
                }

                pap=sum(d,n);
                pap2=n*(n-1)-pap+(n*(pop-n));
                sumll=pap*log(1-exp(-beta))-beta*(pap2);
            }
        }
        if(counter>0 && counter%samgap==0) 
        {   run++;
            R0=R0+pop*beta;
            mean_beta=mean_beta+beta;
            printf(".");
            fprintf(fd,"%f %f %f %f\n",denom[0],sumll,beta,pop*beta);
        }
    }
    printf("\nbeta mean=%f,R0=pop*beta mean=%f\n",(mean_beta/run),(R0/run));
    printf("from  the %d we accept the %d , α=%f\n",(-burn+looplen),miko,((double)(((double)(miko))/((double)(-burn+looplen)))));
    printf("add : A=%d->A2=%d~%f\ndelete B=%d->B2=%d~%f\n",A,A2,((double)(((double)(A2))/((double)(A)))),B,B2,((double)(((double)(B2))/((double)(B)))));
    /*closing the file*/
    fclose(fd);

    /*releasing memory used by the program*/    
    for (i=0; i<n; i++)  
        free(c[i]);
    free(c);             
    
    for (i=0; i<n; i++)  
        free(propc[i]);
    free(propc);            
    
    free(d);            
    free(propd);
    
    return 0;
}

/* the end of program start of function */

int randomIntegern(int n)
{
    int k;
    double number;
    
    /* Generate a random integer number sto(1,n) */
    number = (double)rand()/((double) RAND_MAX+1);
    k=(int)(number*n);
    return (k+1);
}


/* epistrefei apo mia omoiomorfh se [0,1] */

double randomInteger01(void)
{
    double number;
    
    number = (double)rand()/((double) RAND_MAX+1);
    return (number);
}


/* faction pou dinei mia timh apo kanonikh */
double randn (double mu, double sigma)
{
    double U1, U2, W, mult;
    static double X1, X2;
    static int call = 0;
    
    if (call == 1)
    {
        call = !call;
        return (mu + sigma * (double) X2);
    }
    
    do
    {
        U1 = -1 + ((double) rand () / RAND_MAX) * 2;
        U2 = -1 + ((double) rand () / RAND_MAX) * 2;
        W = pow (U1, 2) + pow (U2, 2);
    }
    while (W >= 1 || W == 0);
    
    mult = sqrt ((-2 * log (W)) / W);
    X1 = U1 * mult;
    X2 = U2 * mult;
    
    call = !call;
    
    return (mu + sigma * (double) X1);
}

void loglikelihood(double * result, double beta,int **c,int *d,int n,int pop)
{
    double llike1,sum1,sum2;
    int llike2,i,j;
    llike1=0;    /*arxikopoihsh */
    llike2=0;
    
    if(check(c,d,1,n)==1)    /*an to grafhma einai dekto tote ksekiname*/
    {
        llike2=1;                /* it's valid */
        sum1=0;
        for(i=0;i<n; i++)       /* for all individual */
        {
            sum2=0;
            for(j=0; j<=d[i]; j++)  /* for all links */
            {
                sum2=sum2+choose(d[i],j)*pow((-1),(d[i]-j))*exp(-beta*(pop-j));
            }
            if(sum2<=0)         
            {
                llike2=0;  
            }
            else
            {
                sum1=sum1+log(sum2);  
            }
        }
        llike1=sum1;
    }
    
    result[0]=llike1;   /*return the value of the likelihood if calculated*/
    result[1]=llike2;   /*return 0 if the graph is not valid, 1 with likelihood in result[0] else*/
}

int check(int **c,int *d,int start,int n)
{
    int checkout,i,*hit;
    
    hit=(int *) malloc(n*sizeof(int));
    
    for(i=0; i<n; i++)   /*arxikopoihsh tou dianusmatos hit, einai ena dianusma kai deixnei poios exei stipoithei apo thn asthenia*/
    {
        hit[i]=0;
    }
    hit[0]=1;
    dohits(c,d,start,n,hit);
    if(sum(hit,n)==n)
    {
        checkout=1;
    }
    else
    {
        checkout=0;
    }
    
    free(hit);
    
    return checkout;
}


void dohits(int **c,int *d, int start,int n,int *hit)
{
    int i; 
    if(d[start-1]>0)
    {
        for(i=0; i<d[start-1]; i++)  /*gia to kathe link */
        {
            if(hit[c[start-1][i]-1]==0)  /*follow the link only if it is its first hit */
            {
                hit[c[start-1][i]-1]=1;
                dohits(c,d,c[start-1][i],n,hit); /*here we refer to the actual name of the node we point to, not the position in c so no -1*/
            }
        }
    }
}

long int choose(int n,int k)
{
    double result;
    int i;
    
    if(k>n) return(0);
    
    if(k>n/2) k=n-k;
    
    result=1;
    for(i=1;i<=k;i++)
        result*=(double)(n+1-i)/(i);
    
    return ((long int) result);
}
/* athrisma enos pinaka array */

int sum(int *pa, int num_elements)
{
    int i, sum1=0;
    for (i=0; i<num_elements; i++)
    {
        sum1 = sum1 + pa[i];
    }
    return(sum1);
}

void display_mat_deg(int ** c, int * d, int n)
{
    int i,j;
    
    for(i=0; i<n; i++)       /* display c*/
    {
        printf("(%d)->|",i+1);
        for(j=0; j<n; j++)
            printf("%d ",c[i][j]);
        printf("<->%d",d[i]);
        printf("\n");
   	}
}


/*Function to evaluate mgf = E[exp(-beta I)] 

void mgf(double result, double b )
{
   mgout =exp(-b);
}                         */
 
void award(int*A,int**B,int*A2,int**B2,int n)
{
    int i,j;
    for(i=0; i<n; i++)           
	{
	     A[i]=A2[i];
      	     for(j=0; j<n; j++)
             {
                  B[i][j]=B2[i][j];  
             }
        }
}

int look2(int*vector,int**matrix,int n,int j)
{
     int i;
     printf("\nthe d vector is -->:\n");
     for(i=0; i<n;i++)
     {  
          printf(" %d ",vector[i]);
     }
     printf("\nand the %d's matrix row is:\n",j);       
     for(i=0; i<n; i++)
          {
              printf(" %d ",matrix[j-1][i]);
          }
          printf("\n");   

return 0;
}
