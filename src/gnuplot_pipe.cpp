#include "gnuplot_pipe.h"

FILE *gpc_init_image ()
{
  FILE *pipe;
	
  pipe = popen ("gnuplot > /dev/null 2>&1", "w");                  // Open pipe to Gnuplot and check for error
  if (pipe == NULL)
  {
    printf ("\nCan not find the required Gnuplot executable.");
    printf ("\nGraph creation failure\n");
    exit (1);
  }
  
  fflush (pipe);                                    // flush the pipe
  return (pipe);
}

int gpc_plot_image (FILE *pipe, vector<vec> vectors_A, vec w, int indexToPlot, int bosonOffset = 0)
{
  int i, j;
  char title[80];
  fprintf (pipe, "set term %s size 1200,800 noraise \n", _TERMINAL); // Set the plot
  
  sprintf(title,"A(omega)");
  
  fprintf (pipe, "plot '-' u 1:2 w lp\n"); // Set plot format

  printf("alpha: %d\n", indexToPlot);
	for(int ii=bosonOffset; ii<vectors_A[indexToPlot].size(); ii++) {
    fprintf(pipe,"%5.9f %4.8f\n", w[ii], vectors_A[indexToPlot][ii-bosonOffset]);
    //printf("%5.9f %4.8f\n", w[ii], vectors_A[indexToPlot][ii-bosonOffset]);
  }
  //exit(1);     
	fprintf(pipe, "e\n");                       // End of spectrogram dataset
  fflush(pipe);                                    // Flush the pipe
  //fflush(stdout);                                    // Flush the pipe
  return (0);
}

void gpc_close (FILE *pipe)
{
  fprintf (pipe, "exit\n");                         // Close GNUPlot
  pclose (pipe);                                    // Close the pipe to Gnuplot
}

