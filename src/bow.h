// Boinc tiny interface

#ifdef BOINC

void bow_init(void);
int  bow_get_checkpoint_seconds(int defval);
void bow_report_progress(double fraction_done);
void bow_app_checkpointed(double fraction_done);
int  bow_standalone(void);
void bow_send_trickle_up(const char *variety, double ratio_done);
void bow_finish(int code);
unsigned bow_poll_events(void);

#define BOW_EVENT_QUIT_NORMAL   0x01
#define BOW_EVENT_QUIT_HBT      0x02
#define BOW_EVENT_ABORT         0x04
#define BOW_EVENT_SUSPENDED     0x08

#else

#define bow_init()
#define bow_get_checkpoint_seconds(defval) (defval)
#define bow_finish(code) exit(code)

#endif
