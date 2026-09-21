%module effsource_legacy_elliptic

/* As effsource_legacy, but linking legacy_elliptic/kerr-equatorial.c, in which
   K and E are evaluated from the complementary parameter by AGM. Separate
   shared object, so its global particle state is independent of
   effsource_legacy's and both may be live at once. */

%include "effsource_legacy_api.i"
