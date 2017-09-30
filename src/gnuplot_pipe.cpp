#include "gnuplot_pipe.h"

#define XMIN -50.
#define XMAX 50.
#define YMIN 0.
#define YMAX 1.8


FILE *gpc_init_image ()
{
  FILE *pipe;
	
  pipe = popen ("gnuplot > /dev/null 2>&1", "w");                  // Open pipe to Gnuplot and check for error
  fflush(pipe);                                   // flush the pipe
  if (pipe == NULL)
  {
    printf ("\nCan not find the required Gnuplot executable.");
    printf ("\nGraph creation failure\n");
    exit (1);
  }
  fprintf(pipe,"set xrange[%f:%f]\nset yrange [%f:%f]\n",XMIN,XMAX,YMIN,YMAX);
  fprintf(pipe,"set term %s size 1200,600 noraise \n", _TERMINAL); // Set the plot
  fprintf(pipe,"set tmargin at screen 0.1\nset lmargin at screen 0.1\nset rmargin at screen 0.95\nset bmargin at screen 0.95\n");
  
  fprintf(pipe,"set object 1 rect from %f,%f to %f,%f lw 1 fs empty border lc rgb '#000000' \n", XMIN+0.7*(XMAX-XMIN),YMIN+0.6*(YMAX-YMIN),XMAX,YMAX);
  fprintf(pipe,"plot 0 w l t ''\n"); // Set plot format
  fprintf(pipe,"set label 1 'log {/Symbol c}^2' at %f,%f rotate by 90\n", XMIN+0.68*(XMAX-XMIN), YMIN+0.82*(YMAX-YMIN));
  fprintf(pipe,"set label 2 'log {/Symbol a}' at %f,%f \n",               XMIN+0.83*(XMAX-XMIN), YMIN+0.58*(YMAX-YMIN));
  fprintf(pipe,"set xlabel '{/Symbol w}' \n");
  fprintf(pipe,"set ylabel 'A({/Symbol w})' \n");
  fflush(pipe);                                    // flush the pipe
  return (pipe);
}

int gpc_plot_image (FILE *pipe, vector<vec> vectors_A, vec w, vec lalpha_vec, vec lchi2_vec, double lalpha_min, int indexToPlot, int bosonOffset = 0)
{
  int i, j;
  char title[80];
  
  int sz1= vectors_A[indexToPlot].size();
  int sz2= lalpha_vec.size();
  
  //fprintf(pipe,"set multiplot\nunset logscale\n");
  fprintf (pipe, "plot '-' u 1:2 w lp t '', '-' u 3:4 w p t ''\n");//, '-' u 1:3 w l, '-' u 1:4 w l, '-' u 1:5 w l\n"); // Set plot format
  printf("alpha: %d   %d %d\n", indexToPlot, sz1, sz2);
  
  double x1,y1;
  double x2,y2;
  for(int ii=0; ii<sz2; ii++) {
    if(ii < sz1){
      x1=w[ii];
      y1=vectors_A[indexToPlot][ii];
    }
    else{
      x1=10000.0; // dummy false points
      y1=0.0;
    }
    if(ii < sz2){
      x2=(lalpha_vec[ii]-lalpha_min)/(lalpha_vec[0]-lalpha_min)*0.3*(XMAX-XMIN)+XMIN+0.7*(XMAX-XMIN);       // complicated tricks to plot the curve on the same plot, but in the corner.
      y2=lchi2_vec[ii]/lchi2_vec[0]*0.4*(YMAX-YMIN)+YMIN+0.55*(YMAX-YMIN);
    }
    else{
      x2=1000; // dummy false points
      y2=1000;
    }
    fprintf(pipe,"%1.3f %1.3f %1.3f %1.3f\n", x1,y1,x2,y2);
  }
  
      

    //fprintf(pipe,"%5.9f %4.8f %4.8f %4.8f %4.8f\n", w[ii], vectors_A[indexToPlot][ii-bosonOffset], vectors_A[indexToPlot][ii-bosonOffset]*1.1, vectors_A[indexToPlot][ii-bosonOffset]*1.3, vectors_A[indexToPlot][ii-bosonOffset]*1.2);
    //printf("%5.9f %4.8f\n", w[ii], vectors_A[indexToPlot][ii-bosonOffset]);
  
  //exit(1);     
  fprintf(pipe, "e\n\n");                       // End of spectrogram dataset
  
  //fprintf(pipe,"unset yrange\nset format y ''\nunset ylabel\nset format x ''\nunset xlabel\nunset border\n");
  //fprintf(pipe,"set yrange [2.0:%1.2f]\nset xrange [%1.2f:%1.2f]\n", lchi2_vec(0), lalpha_min, lalpha_vec(0));
  //fprintf(pipe,"set tmargin at screen 0.6\nset lmargin at screen 0.7\nset rmargin at screen 0.95\nset bmargin at screen 0.95\n");
  //fprintf (pipe, "plot '-' u 1:2 w l \n"); // Set plot format
  //printf("alpha: %d\n", indexToPlot);
  //for(int ii=0; ii<lalpha_vec.size(); ii=ii+5) {
  //  fprintf(pipe,"%5.9f %4.8f\n", lalpha_vec[ii]+20.0, lchi2_vec[ii]/lchi2_vec[0] );
    //printf("%5.9f %4.8f\n", w[ii], vectors_A[indexToPlot][ii-bosonOffset]);
  //}
  
  //fprintf(pipe, "e\n\n\n");                       // End of spectrogram dataset
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

