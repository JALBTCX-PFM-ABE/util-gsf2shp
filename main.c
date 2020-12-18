
/*********************************************************************************************

    This is public domain software that was developed by or for the U.S. Naval Oceanographic
    Office and/or the U.S. Army Corps of Engineers.

    This is a work of the U.S. Government. In accordance with 17 USC 105, copyright protection
    is not available for any work of the U.S. Government.

    Neither the United States Government, nor any employees of the United States Government,
    nor the author, makes any warranty, express or implied, without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE, or assumes any liability or
    responsibility for the accuracy, completeness, or usefulness of any information,
    apparatus, product, or process disclosed, or represents that its use would not infringe
    privately-owned rights. Reference herein to any specific commercial products, process,
    or service by trade name, trademark, manufacturer, or otherwise, does not necessarily
    constitute or imply its endorsement, recommendation, or favoring by the United States
    Government. The views and opinions of authors expressed herein do not necessarily state
    or reflect those of the United States Government, and shall not be used for advertising
    or product endorsement purposes.

*********************************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>

#include "shapefil.h"

#include "nvutility.h"

#include "gsf.h"

#include "version.h"


extern int32_t gsfError;


int32_t main (int32_t argc, char **argv)
{
  gsfDataID           id;
  gsfRecords          gsfrecord;
  int32_t             i, j, recnum, percent = 0, old_percent = -1, gsf_hnd, count = 0, tide_type;
  char                gsf_name[512], shape_name[512], prj_file[512];
  SHPHandle           shp_hnd;
  SHPObject           *shape;
  DBFHandle           dbf_hnd;
  FILE                *prj_fp;
  double              *x_array = NULL, *y_array = NULL, *z_array = NULL, *m_array = NULL, lat[3], lon[3], ang1, ang2, lateral, 
                      ping_time, dep;


  printf("\n\n %s \n\n",VERSION);


  if (argc < 3)
    {
      fprintf (stderr, "This program converts SWATH_BATHYMETRY_PING records in a GSF file\n");
      fprintf (stderr, "to MultiPointZ records in an ESRI shape file.  The Z value is depth in\n");
      fprintf (stderr, "meters.  The M (measure) value is time in seconds from 00:00:00 on\n");
      fprintf (stderr, "January 1, 1970, Coordinated Universal Time (UTC) (POSIX.1 time).\n");
      fprintf (stderr, "SHAPE_FILE_BASE_NAME must not exceed 8 characters.\n");
      fprintf (stderr, "Four files will be created with .shp, .shx, .dbf, and .prj extensions.\n\n");
      fprintf (stderr, "USAGE: gsf_shape GSF_FILE SHAPE_FILE_BASE_NAME\n\n");
      exit (-1);
    }


  strcpy (gsf_name, argv[1]);

  if (gsfOpen (gsf_name, GSF_READONLY_INDEX, &gsf_hnd))
    {
      gsfPrintError (stderr);
      exit (-1);
    }

  printf("File : %s\n\n", gsf_name);


  strcpy (shape_name, argv[2]);

  if ((shp_hnd = SHPCreate (shape_name, SHPT_MULTIPOINTZ)) == NULL)
    {
      perror (shape_name);
      exit (-1);
    }


  if ((dbf_hnd = DBFCreate (shape_name)) == NULL)
    {
      perror (shape_name);
      exit (-1);
    }

  if (DBFAddField (dbf_hnd, "t_type", FTString, 10, 0) == -1)
    {
      perror ("Adding field to dbf file");
      exit (-1);
    }


  if (DBFAddField (dbf_hnd, "t_corr", FTDouble, 6, 2) == -1)
    {
      perror ("Adding field to dbf file");
      exit (-1);
    }


  //  Stupid freaking .prj file

  strcpy (prj_file, shape_name);
  strcat (prj_file, ".prj");

  if ((prj_fp = fopen (prj_file, "w")) == NULL)
    {
      perror (prj_file);

      exit (-1);
    }

  fprintf (prj_fp, "GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563,AUTHORITY[\"EPSG\",\"7030\"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY[\"EPSG\",\"6326\"]],PRIMEM[\"Greenwich\",0,AUTHORITY[\"EPSG\",\"8901\"]],UNIT[\"degree\",0.01745329251994328,AUTHORITY[\"EPSG\",\"9108\"]],AXIS[\"Lat\",NORTH],AXIS[\"Long\",EAST],AUTHORITY[\"EPSG\",\"4326\"]]\n");
  fclose (prj_fp);


  recnum = 0;
  j = 0;
  while (1)
    {
      id.recordID = GSF_RECORD_SWATH_BATHYMETRY_PING;
      id.record_number = recnum;
      if (gsfRead (gsf_hnd, GSF_RECORD_SWATH_BATHYMETRY_PING, &id, &gsfrecord, NULL, 0) < 0)
        {
          if (gsfError == GSF_INVALID_RECORD_NUMBER) break;
        }


      lat[0] = gsfrecord.mb_ping.latitude;
      lon[0] = gsfrecord.mb_ping.longitude;
      ping_time = (double) gsfrecord.mb_ping.ping_time.tv_sec + (double) gsfrecord.mb_ping.ping_time.tv_nsec *
        1000000000.0;


      if ((x_array = (double *) calloc (gsfrecord.mb_ping.number_beams, sizeof (double))) == NULL)
        {
          perror ("Allocating X array");
          exit (-1);
        }

      if ((y_array = (double *) calloc (gsfrecord.mb_ping.number_beams, sizeof (double))) == NULL)
        {
          perror ("Allocating Y array");
          exit (-1);
        }

      if ((z_array = (double *) calloc (gsfrecord.mb_ping.number_beams, sizeof (double))) == NULL)
        {
          perror ("Allocating Z array");
          exit (-1);
        }

      if ((m_array = (double *) calloc (gsfrecord.mb_ping.number_beams, sizeof (double))) == NULL)
        {
          perror ("Allocating M array");
          exit (-1);
        }


      ang1 = gsfrecord.mb_ping.heading + 90.0;
      ang2 = gsfrecord.mb_ping.heading;


      j = 0;


      /*
       * Make sure position and ping is valid.   
       */

      if ((lat[0] <= 90.0) && (lon[0] <= 180.0) && !(gsfrecord.mb_ping.ping_flags & GSF_IGNORE_PING))
        {
          for (i = 0 ; i < gsfrecord.mb_ping.number_beams ; i++) 
            {
              if (gsfrecord.mb_ping.depth == NULL)
                {
                  dep = gsfrecord.mb_ping.nominal_depth[i];
                }
              else
                {
                  dep = gsfrecord.mb_ping.depth[i];
                }

              if (dep != 0.0 && gsfrecord.mb_ping.beam_flags != NULL && 
                  !(check_flag (gsfrecord.mb_ping.beam_flags[i], NV_GSF_IGNORE_NULL_BEAM))) 
                {
                  /*
                   * Adjust for cross track position.  
                   */

                  lateral = gsfrecord.mb_ping.across_track[i];
                  newgp (lat[0], lon[0], ang1, lateral, &lat[1], &lon[1]);


                  /*
                   * if the along track array is present then use it 
                   */
                
                  if (gsfrecord.mb_ping.along_track != (double *) NULL) 
                    {
                      lateral = gsfrecord.mb_ping.along_track[i];

                      newgp (lat[1], lon[1], ang2, lateral, &lat[2], &lon[2]);

                      x_array[j] = lon[2];
                      y_array[j] = lat[2];
                    }
                  else
                    {
                      x_array[j] = lon[1];
                      y_array[j] = lat[1];
                    }
                  z_array[j] = dep;
                  m_array[j] = ping_time;

                  j++;
                }
            }
        }

      if (j)
        {
          shape = SHPCreateObject (SHPT_MULTIPOINTZ, -1, 0, NULL, NULL, j, x_array, y_array, z_array, m_array);

          SHPWriteObject (shp_hnd, -1, shape);
          SHPDestroyObject (shape);

          tide_type = (gsfrecord.mb_ping.ping_flags & (GSF_PING_USER_FLAG_15 | GSF_PING_USER_FLAG_14)) >> 14;

          switch (tide_type)
            {
            case 0:
            default:
              DBFWriteStringAttribute (dbf_hnd, count, 0, "None");
              break;

            case 1:
              DBFWriteStringAttribute (dbf_hnd, count, 0, "Predicted");
              break;

            case 2:
              DBFWriteStringAttribute (dbf_hnd, count, 0, "Observed");
              break;

            case 3:
              DBFWriteStringAttribute (dbf_hnd, count, 0, "Verified");
              break;
            }

          DBFWriteDoubleAttribute (dbf_hnd, count, 1, (double) gsfrecord.mb_ping.tide_corrector);
          count++;
        }

      free (x_array);
      free (y_array);
      free (z_array);
      free (m_array);


      recnum++;


      percent = gsfPercent(gsf_hnd);
      if ( old_percent != percent )
        {
          printf ("%3d%% processed    \r", percent);
          fflush (stdout);
          old_percent = percent;
        }
    }


  percent = 100;
  printf ("%3d%% processed    \n", percent);
        
  gsfClose(gsf_hnd);
  printf("\n");


  SHPClose (shp_hnd);
  DBFClose (dbf_hnd);


  return (0);
}
