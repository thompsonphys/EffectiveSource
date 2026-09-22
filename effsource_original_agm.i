%module effsource_original_agm

/* As effsource_original, but linking original_agm/kerr-equatorial.c, in which
   K and E are evaluated from the complementary parameter by AGM. Separate
   shared object, so its global particle state is independent of
   effsource_original's and both may be live at once. */

%include "effsource_original_api.i"
