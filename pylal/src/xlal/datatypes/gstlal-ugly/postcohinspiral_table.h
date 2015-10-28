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
  LIGOTimeGPS	end_time;
  INT4		is_background;
  INT4		livetime;
  CHAR		ifos[MAX_ALLIFO_LEN];
  CHAR		pivotal_ifo[MAX_IFO_LEN];
  INT4		tmplt_idx;
  INT4		pix_idx;
  REAL4		maxsnglsnr;	
  REAL4         cohsnr;
  REAL4         nullsnr;
  REAL4         chisq;
  REAL4		spearman_pval;
  REAL4		fap;
  REAL4		far;
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
  REAL8		ra;
  REAL8		dec;
}
PostcohInspiralTable;
#endif /* __POSTCOH_TABLE_H */
