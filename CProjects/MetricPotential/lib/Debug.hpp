#ifndef H_DEBUG_H
#define H_DEBUG_H

int countCalls(int);

#ifdef DEBUG
#include <cstdio>
#define ECHOS printf("[%s] START\n", __PRETTY_FUNCTION__)
#define ECHOF printf("[%s] END\n", __PRETTY_FUNCTION__)
#else
#define ECHOS
#define ECHOF
#endif

#ifdef COUNT
#include <cstdio>
#define COUNT_B(ch)  countCalls(ch)
#define COUNT_P(ch) printf("[%s] ch:%d, cn:%d\n",__PRETTY_FUNCTION__, ch, countCalls(ch));
#else
#define COUNT_B
#define COUNT_P
#endif

#endif
