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

int gpc_plot_image (FILE *pipe, vector<vec> vectors_A, vec w, vec alpha_vec, vec chi2_vec, double lalpha_min, int indexToPlot, int bosonOffset = 0)
{
  int i, j;
  char title[80];
  fprintf (pipe, "set term %s size 1200,600 noraise \n", _TERMINAL); // Set the plot
  
  //fprintf(pipe,"set multiplot\nunset logscale\n");
  fprintf(pipe,"set yrange [0:1.8]\nset xrange[-50:50]\n");
  fprintf(pipe,"set tmargin at screen 0.1\nset lmargin at screen 0.1\nset rmargin at screen 0.95\nset bmargin at screen 0.95\n");
  fprintf (pipe, "plot '-' u 1:2 w lp\n"); // Set plot format
  printf("alpha: %d\n", indexToPlot);
  for(int ii=bosonOffset; ii<vectors_A[indexToPlot].size(); ii++) {
    fprintf(pipe,"%5.9f %4.8f\n", w[ii], vectors_A[indexToPlot][ii-bosonOffset]);
    //printf("%5.9f %4.8f\n", w[ii], vectors_A[indexToPlot][ii-bosonOffset]);
  }
  //exit(1);     
  fprintf(pipe, "e\n\n\n");                       // End of spectrogram dataset
  
  //fprintf(pipe,"unset yrange\nset format y ''\nunset ylabel\nset format x ''\nunset xlabel\nunset border\n");
  //fprintf(pipe,"set yrange [2.0:%1.2f]\nset xrange [%1.2f:%1.2f]\n", chi2_vec(0), lalpha_min, alpha_vec(0));
  //fprintf(pipe,"set tmargin at screen 0.6\nset lmargin at screen 0.7\nset rmargin at screen 0.95\nset bmargin at screen 0.95\n");
  fprintf (pipe, "plot '-' u 1:2 w l \n"); // Set plot format
  //printf("alpha: %d\n", indexToPlot);
  for(int ii=0; ii<alpha_vec.size(); ii=ii+5) {
    fprintf(pipe,"%5.9f %4.8f\n", alpha_vec[ii]+20.0, chi2_vec[ii]/chi2_vec[0] );
    //printf("%5.9f %4.8f\n", w[ii], vectors_A[indexToPlot][ii-bosonOffset]);
  }
  
  fprintf(pipe, "e\n\n\n");                       // End of spectrogram dataset
  //fprintf(pipe, "unset multiplot\n");
  
  fflush(pipe);      
  //std::cin.ignore();                              // Flush the pipe
  return (0);
}

void gpc_close (FILE *pipe)
{
  fprintf (pipe, "exit\n");                         // Close GNUPlot
  pclose (pipe);                                    // Close the pipe to Gnuplot
}

