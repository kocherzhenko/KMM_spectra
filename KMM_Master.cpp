// Master file: program for running all linear KMM calculations - (c) Aleksey Kocherzhenko, 2014

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include "Path.h"

FILE *fp;

// Start of the main program
int main (int argc, char **argv)
{
  int i, j;
  int units, fourthorder;
  double coup;
  char c;
  char order[256];
  char *location, *backlocation, *buff;
  char *command;

  printf("Starting KMM linear spectroscopy calculation...\n\n");
  
  if (argc == 1)
  {
    printf("Error: not enough arguments (please enter input file name)\nExiting...\n\n");
    return 1;
  }
  else if (argc > 2)
  {
    printf("Error: too many arguments\nExiting...\n\n");
    return 1;
  }
  else
  {
    // read # of sites
    fp = fopen(argv[1], "r");
    if (fp == NULL)
    {
      printf("Input file %s does not exist. Exiting...\n", argv[1]);
      return 1;
    }
    do { c = fgetc(fp); } while (c != ':');
    j = fscanf(fp, " %d ", &units);
    if (j < 1)
    {
      printf("Error reading input file %s, exiting...\n", argv[1]);
      return 1;
    }
    // read coupling (in units of epsilon)
    do { c = fgetc(fp); } while (c != ':');
    j = fscanf(fp, " %lf ", &coup);
    if (j < 1)
    {
      printf("Error reading input file %s, exiting...\n", argv[1]);
      return 1;
    }
    do { c = fgetc(fp); } while (c != ':');
    j = fscanf(fp, " %s ", order);
    if (j < 1)
    {
      printf("Error reading input file %s, exiting...\n", argv[1]);
      return 1;
    }
    else
    {
      if ((*order == 'y') || (*order == 'Y')) { fourthorder = 1; }
      else if ((*order == 'n') || (*order == 'N')) { fourthorder = 0; }
      else
      {
        printf("Error reading input file %s, exiting...\n", argv[1]);
        return 1; 
      }
    }
    fclose(fp);
  }
  
  if (units % 2 != 0)
  {
    printf("Number of sites should be even, exiting...\n");
    return 1;
  }

  if ((fabs(coup)>0.499999) && (fabs(coup)<0.500001))
  {
    printf("Coupling value too close to 0.5: divergent results, exiting...\n");
    return 1;
  }

  if (coup > 0.0)
  {
    printf("Coupling value should be negative, exiting...\n");
    return 1;
  }
  
//  printf(" %d  %10.5f  %d  %d\n", units, coup, fourthorder, *order);

// Go to directory where the excecutables are
  backlocation = getcwd(buff, 256);
//  printf("%s \n", backlocation);
  location = getenv("HOME");
  location = strcat(location, LOCATION);
  i=chdir(location);
//  printf(" %d %d \n", i, errno);
  
// Write "Input.in"
  fp = fopen("Input.in", "w");
  fprintf(fp, "%6d\n", units);
  fprintf(fp, "%+15.10f\n", coup);
  fprintf(fp, "%6d\n", fourthorder);
  fclose(fp);
  
//  system("sh srccomp");
  system("./matrixSetup");
  system("./finiteEqs");
  system("./vacCalc");
  system("./weakParams");
  system("./lineInt");
  system("./higherStrong");
  system("./minMax");
  system("./gridConvert");
  
  system("mv *.out scratch/text"); 
  system("mv *.in scratch/text"); 
  
  printf("Copying output files from scratch directory...\n");
  i = chdir("scratch/text");
//  location = getcwd(buff, 256);
//  printf(" %s \n", location);
  
  command = (char *) malloc((strlen(backlocation)+5)*sizeof(char));
  command[0] = 'm';
  command[1] = 'v';
  command[2] = ' ';
  command[3] = '*';
  command[4] = ' ';
  command = strcat(command, backlocation);
//  printf(" %s \n", command);
  system(command);
  i = chdir(backlocation);
//  location = getcwd(buff, 256);
//  printf(" %s \n", location);  
  free(command);
  
  system("rm *.in");
  system("rm FinKsVsAlpha.out");
  system("rm InfKsVsAlpha.out");
  system("rm XEigvals.out");
 
  if (fabs(coup) < 0.5)
  {
    system("rm 1InfIntens.out");
    system("rm 3InfIntens.out");
    system("rm FinK123s.out");
    system("rm InfK123s.out");
    system("rm FinJsVsAlpha.out");
    system("rm InfJsVsAlpha.out");
    system("rm Overlap.out");
  }
  else
  {
    system("rm 2InfIntens.out");
    if (fourthorder == 1) { system("rm FinK1234s.out"); }
  }
  
  printf("KMM linear spectroscopy calculation finished!\n\n");

  return 0;
}

