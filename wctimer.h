/* ******************************************************/
/*                                                      */
/*  timer.h - timing structs and inlines                */
/*                                                      */
/*  created: 6/30/19                                    */
/*                                                      */
/* ******************************************************/

#include <inttypes.h>
#include <time.h>
#include <unistd.h>

// compiler intrinsics for TSC timing, should be okay on clang, gcc, and icc
#include <x86intrin.h>

/* primary timing data structure */
struct wc_timer_s {
  uint64_t total;  // total time elapesd
  uint64_t last;   // timestamp when timer was started
  uint64_t temp;   // used to accumulate times
};
typedef struct wc_timer_s wc_timer_t;


/*
 * wc_get_wtime() -- use system utilities to get timing info
 */
static inline struct timespec wc_get_wtime() {
    struct timespec ts;
#if __APPLE__
      clock_gettime(_CLOCK_MONOTONIC, &ts);
#else
      clock_gettime(CLOCK_MONOTONIC_RAW, &ts);
#endif // clock_gettime
          return ts;
}


/*
 * wc_get_tsctime() -- use assembly instruction to get TSC counter value
 */
static inline uint64_t wc_get_tsctime() {
  return __rdtsc();
}


/*
 * macros for init/start/stop and reading timer values
 */
#define WC_CPU_HZ                 tsc_cpu_hz*(double)1e6
#define WC_INIT_TIMER(TMR)       do { TMR.total = 0; TMR.temp = 0; TMR.last = 0; } while (0)
#define WC_START_TIMER(TMR)      TMR.last = wc_get_tsctime();
#define WC_STOP_TIMER(TMR)       do {\
                                      TMR.temp = wc_get_tsctime();\
                                      TMR.total += TMR.temp - TMR.last;\
                                    } while (0)

#define WC_READ_TIMER(TMR)          TMR.total
#define WC_READ_TIMER_NSEC(TMR)     (double)(WC_READ_TIMER(TMR)/(WC_CPU_HZ)) * (double)1e9
#define WC_READ_TIMER_USEC(TMR)     (double)(WC_READ_TIMER(TMR)/(WC_CPU_HZ)) * (double)1e6
#define WC_READ_TIMER_MSEC(TMR)     (double)(WC_READ_TIMER(TMR)/(WC_CPU_HZ)) * (double)1000.0
#define WC_READ_TIMER_SEC(TMR)      (double)WC_READ_TIMER(TMR)/(WC_CPU_HZ)


extern double tsc_cpu_hz;  // global variable, shhh.

// timer.c prototypes
void wc_tsc_calibrate(void);
