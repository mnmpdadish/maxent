

#include "maxEnt_data.h"

int main(int arg_N, char *args[])
{
	OmegaMaxEnt_data maxent1(arg_N, args);
	maxent1.loop_run();
		
    return 0;
}
