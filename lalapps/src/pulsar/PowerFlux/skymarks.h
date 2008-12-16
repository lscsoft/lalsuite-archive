#ifndef __SKYMARKS_H__
#define __SKYMARKS_H__


typedef struct {
	int type;

	#define SKYMARK_NOP  0
	#define SKYMARK_CLOSEST  1
	#define SKYMARK_DISK  2
	#define SKYMARK_BAND  3
	#define SKYMARK_RESPONSE  4
	#define SKYMARK_ECHO_RESPONSE  5
	#define SKYMARK_LINE_RESPONSE  6

	#define SKYMARK_END 255

	#define SKYMARK_ERROR -1

	int band_to;
	int band_from;


	union {
		struct {
			float ra;
			float dec;
			/* precomputed values */
			int closest_point;
			} closest;

		struct {
			float ra;
			float dec;
			float radius;
			/* precomputed values */
			float cos_radius;
			int closest_point;
			} disk;
		struct {
			float ra;
			float dec;
			float level1;
			float level2;
			/* precomputed values */
			float x0;
			float y0;
			float z0;
			} band;
		struct {
			float ra;
			float dec;
			float weight_ratio_level;
			float bin_tolerance;
			float spindown_tolerance;
			} response;
		struct {
			float ra;
			float dec;
			float ref_spindown;
			float weight_ratio_level;
			float bin_tolerance;
			float spindown_tolerance;
			} echo_response;
		struct {
			float weight_ratio_level;
			float bin_tolerance;
			} line_response;
		
		} p;
	} SKYMARK;


SKYMARK * compile_marks(SKY_GRID *sky_grid, char *s, int length);
int mark_sky_point(SKYMARK *sm, int point, float ra, float dec, float spindown1);

#endif
