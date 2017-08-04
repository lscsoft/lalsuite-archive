/* 
 * Copyright (C) 2014 Qi Chu <qi.chu@ligo.org>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

/*
# This should be the place we define our own table: postcoh_inspiral
# Compare it with the tables in lalmetaio/ src/ LIGOMetadataTables.h
*/

#ifndef __POSTCOH_TABLE_H__
#define __POSTCOH_TABLE_H__

#include <lal/LALStdlib.h> // for the datatypes
#include <lal/Date.h> // for the LIGOTimeGPS

#define MAX_IFO_LEN 4 
#define MAX_ALLIFO_LEN 14 
#define MAX_SKYMAP_FNAME_LEN 50

typedef struct
tagPostcohInspiralTable
{
  struct tagPostcohInspiralTable *next;
  long		process_id;
  long		event_id;
  LIGOTimeGPS	end_time;
  LIGOTimeGPS	end_time_L;
  LIGOTimeGPS	end_time_H;
  LIGOTimeGPS	end_time_V;
  INT4		is_background;
  INT4		livetime;
  CHAR		ifos[MAX_ALLIFO_LEN];
  CHAR		pivotal_ifo[MAX_IFO_LEN];
  INT4		tmplt_idx;
  INT4		pix_idx;
  REAL4		snglsnr_L;
  REAL4		snglsnr_H;
  REAL4		snglsnr_V;
  REAL4		coaphase_L;
  REAL4		coaphase_H;
  REAL4		coaphase_V;
  REAL4		chisq_L;
  REAL4		chisq_H;
  REAL4		chisq_V;
  REAL4         cohsnr;
  REAL4         nullsnr;
  REAL4         cmbchisq;
  REAL4		spearman_pval;
  REAL4		fap;
  REAL4		far_h;
  REAL4		far_l;
  REAL4		far_v;
  REAL4		far_h_1w;
  REAL4		far_l_1w;
  REAL4		far_v_1w;
  REAL4		far_h_1d;
  REAL4		far_l_1d;
  REAL4		far_v_1d;
  REAL4		far_h_2h;
  REAL4		far_l_2h;
  REAL4		far_v_2h;
  REAL4		far;
  REAL4		far_2h;
  REAL4		far_1d;
  REAL4		far_1w;
  CHAR		skymap_fname[MAX_SKYMAP_FNAME_LEN];			// location of skymap
  REAL8		template_duration;
  REAL4		mass1;
  REAL4		mass2;
  REAL4		mchirp;
  REAL4		mtotal;
  REAL4		spin1x;
  REAL4		spin1y;
  REAL4		spin1z;
  REAL4		spin2x;
  REAL4		spin2y;
  REAL4		spin2z;
  REAL4		eta;
  REAL8		ra;
  REAL8		dec;
  REAL8		deff_L;
  REAL8		deff_H;
  REAL8		deff_V;
}


PostcohInspiralTable;
#endif /* __POSTCOH_TABLE_H */
