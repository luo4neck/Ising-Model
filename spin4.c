#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#define JKbT 0.4
int delE(int lts[40][40], int i, int j, int dir)
{
    int up=0, dn=0, rt=0, lt=0, sum=0, rsum=0, rdir=0, d=0;
    if(dir == 1)
        rdir = -1;
    else 
        rdir =  1;
    if(i == 0)
        {up = lts[39][j], dn = lts[i+1][j];}
    else if(i == 39)
        {up = lts[i-1][j],  dn = lts[0][j];}
    else 
        {up = lts[i-1][j],dn = lts[i+1][j];}

    if(j == 0)
        {lt = lts[i][39], rt = lts[i][j+1];}
    else if(j == 39)
        {lt = lts[i][j-1],  rt = lts[i][0];}
    else 
        {lt = lts[i][j-1],rt = lts[i][j+1];}

    sum=dir*up + dir*dn + dir*lt + dir*rt;
    rsum=rdir*up + rdir*dn + rdir*lt + rdir*rt;

    return d=sum-rsum;
}

void print(int lts[40][40])
{
    int i,j;
    for(i=0;i<40;i++)
    {
        for(j=0;j<40;j++)
        {
            if( lts[i][j] > 0)
                printf("+,");
            else
                printf("-,");
        }
        printf("\n");
    }
    printf("\n");
}

int avg(int lts[40][40])
{
    int i,j,sum=0;
    for(i=0;i<40;i++)
    {
        for(j=0;j<40;j++)
            sum=sum + lts[i][j];
    }
    return sum;
}

int main()
{
    FILE *fp;
    FILE* gp;
    FILE* bp;
    #ifdef WIN32 
        gp = _popen("gnuplot -persist","w"); 
        bp = _popen("gnuplot -persist","w"); 
    #else 
        gp = popen("gnuplot -persist","w"); 
        bp = popen("gnuplot -persist","w"); 
    #endif

    int lts[40][40]={}, Ns=0, t=0, ii=0, i, j;//lts[i][j] is the state, 
    double random, mean, var, rsq=0, sum=0, sumsq=0, M[1000]={0}, R[20]={0}, E1=0, E2=0, tao=0, bin[50]={0};
    for (i=0;i<40;i++)// give 1 to all the elements in the lattice..
    {
        for(j=0;j<40;j++)
            lts[i][j]=1;
    }

    srand48(time(NULL));
    for(Ns=1;Ns<=25000;Ns++)// sweep for 5000 times and then begin measurements..
    {
        for (i=0;i<40;i++)
        {
            for(j=0;j<40;j++)
            {
                if( delE(lts,i,j,lts[i][j]) > 0 )//how to deal with 0??
	        {
	            lts[i][j] = -1 * lts[i][j];
	        }
                else
	        {
	            random = drand48();
	            if ( random < exp(JKbT * delE(lts,i,j,lts[i][j])) )
	            lts[i][j] = -1 * lts[i][j];
	        }
            }
        }
  
        if(Ns > 5000 && Ns%20 == 0) // make measurements and record in array..
        {
            M[ii] = avg(lts) / 1600.0;
            bin[t]=bin[t]+M[ii]; // store M in to a bin..
            if((ii+1)%20 == 0 )//check the time to change to next bin..
                t=t+1;//t is the bin number..
            sum = sum + M[ii];
            sumsq = sumsq + M[ii] * M[ii];
            printf("Measurements No.%d (Ns = %d) average = %f\n",ii+1,Ns, M[ii] );    
            ii=ii+1;
        }
    }

    printf("\nJ/Kb T is: %.2f\n",JKbT);  
    printf("After 1000 times of measurements the spin state of the lattice is: \n");
    printf("'+' means 1 and '-' means -1\n");
    print(lts);

    mean = sum/1000.0;
    var = sumsq/1000.0;
    rsq = sumsq/1000.0 - mean*mean;
    printf("Variance is %f, mean = %f.\n\n",sqrt(var), mean);

/*-----------------------begin to deal with plot part-----------------------*/

    if ( (fp=fopen("plot.dat","w+")) == NULL )// open a dat file to plot
    {
        printf("cannot show the result of bin");
        exit(0);
    }
	
    for(t=1;t<=20;t++)//normalized autocorrelation function part..
    {
        E1 = 0, E2 = 0;
        for(i=0;i<1000;i++)
        {
            if(i < 1000-t)	E1 = E1 + M[i]*M[i+t] - mean*M[i];
            if( i >= t+1 )	E2 = E2 + mean*M[i];
        }
  
	R[t-1] = ((E1 - E2)/(1000.0-(double)t) + mean*mean) / rsq;
  
        E1 = 0; // just reuse the E1..
        for(i=1;i<=t;i++)// integrated autocorrelation function part..
        {
            E1 = E1 + (1 - (double)i / (double)t ) * R[i-1];
        }
        tao = 1 + 2 * E1;
        printf("When t=%d, ρ(%d) = %f, τ(%d) = %f\n",t, t, R[t-1], t, tao);
        fprintf(fp, "%d  %f  %f\n", t, R[t-1], tao);//write on dat file to plot
    }
    fclose(fp);  

    if(gp == NULL)  return-1; 
	
    fprintf(gp,"set isosample 20\n"); 
    fprintf(gp,"set contour\n");
    fprintf(gp,"set xlabel 'x axis is t'\n");
    fprintf(gp,"set ylabel 'green line is ρ and red line is τ'\n"); 
    fprintf(gp,"plot 'plot.dat' u 1:2 w l lt 1 lw 1,'plot.dat' u 1:3 w l lt 2 lw 1\n");
    fprintf(gp,"pause -1\n"); 
    fclose(gp);

/*-----------------begin to deal with binning part--------------*/

    if ( (fp=fopen("bin.dat","w+")) == NULL )//open a dat file for histogram..
    {
        printf("cannot show the result of bin");
        exit(0);
    }
    
    for(i=0;i<50;i++)
    {
        fprintf(fp, "%f  %f\n", (i-24)/25.0, bin[i]);
    }
    fclose(fp);

    if(bp == NULL) return-1;
	
    fprintf(bp,"set contour\n");
    fprintf(bp,"set style data histogram\n");
    fprintf(bp,"set style histogram clustered gap 0.5\n");
    fprintf(bp,"set style fill solid 0.4 border\n");
    fprintf(bp,"set xlabel 'x is number of bins'\n"); 
    fprintf(bp,"set ylabel 'J/KbT = %.2f'\n",JKbT); 
    fprintf(bp,"plot [-1:51] 'bin.dat' u 2\n");
    fprintf(bp,"pause -1\n"); 
    fclose(bp);

    return 0;
}
