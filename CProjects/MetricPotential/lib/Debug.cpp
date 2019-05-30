#include "Debug.hpp"

#define MAX_CHANNEL 64

int countCalls(int channel)
{
    static unsigned long int mem[MAX_CHANNEL] = { 0 };

    if(channel>0){
	return (mem[channel]++);
    }
    else if(channel<0){
	return (mem[-channel] = 0);
    }

    return -1;
}
