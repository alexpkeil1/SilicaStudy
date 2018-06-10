*clear the log window and the output window;
DM LOG 'clear;' CONTINUE; DM OUT 'clear;' CONTINUE; 
/**********************************************************************************************************************
* Author: Alex Keil
* Project: Pooled silica exposed worker study
* Released under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html
* this program estimates survival curves using the pseudo-population data created in the gformula_silica.sas program
**********************************************************************************************************************/

%INCLUDE "gformula_silica.sas"; *first get pseudo-cohort;

LIBNAME an  "&base/data/an";
LIBNAME gformula  "&base/data/derived";


* macro to get survival curves for a specific intervention;
%MACRO mksurvdata (ds=nc, invar=agein, outvar=age, suff=, intervention=,breakties=1);
 *PHREG cannot handle negative times, so convert sas dates to number of days from the minimum start time (will be converted back at end);
 %LET mintime = -900000; *cause an obvious problem in the graph if assignment of this does not work;
 TITLE "gformula.cidata_&outvar._&ds&suff";
 DATA ___sur;
  SET &DS;
  %IF &intervention NE %THEN IF intervention=&intervention;;
 
 PROC MEANS DATA = ___SUR NOPRINT;VAR &invar;
  OUTPUT OUT=___SURMIN MIN=minval;
 DATA _null_;SET ___surmin; CALL SYMPUT('mintime', PUT(minval, BEST9.));RUN;
 
 DATA ___sur;
  SET ___sur;
  &invar =  &invar - &mintime;
  &outvar = &outvar - &mintime;
  d_multi = d_lc + 2*d_other;
  *random breaking of ties for mortality outcomes only (introduced 3/11/2017);
  %IF &breakties=1 %THEN %DO;
  CALL STREAMINIT(12323100);
  IF d_multi>0 THEN DO;
   ct=1;
   __tmp = &outvar;
   DO WHILE (__tmp <= &invar OR ct=1) ;
    __tmp = &outvar - (RAND('uniform'))*0.0001;
    ct=ct+1; 
   END;
   &outvar=__tmp;
   DROP ct;
  END;
  %END;
 RUN;
 *causes of death;
 PROC PHREG DATA = ___sur NOPRINT;
  WEIGHT sampwt / NORM;
  MODEL (&invar &outvar)*d_multi(0) =  / TIES=EFRON;
  OUTPUT OUT=surv_multi_&outvar._&ds(KEEP=d_multi &invar &outvar surv_y) SURVIVAL=surv_y / METHOD=PL ;
 RUN;
 PROC SORT DATA = surv_multi_&outvar._&ds(WHERE=(d_multi>0)); BY &outvar;RUN;
 
 DATA cidata_&outvar._&ds&suff;
  SET surv_multi_&outvar._&ds(KEEP=&outvar surv_y d_multi);
    &outvar = &outvar + &mintime;
 RUN;
 PROC SORT DATA = cidata_&outvar._&ds&suff; BY &outvar DESCENDING surv_y;
 DATA gformula.cidata_&outvar._&ds&suff (KEEP = &outvar surv_y  ci_: mark) test;
  SET cidata_&outvar._&ds&suff;
  BY &outvar DESCENDING surv_y;
  RETAIN lastsurv_y 1 ci_ac ci_lc  ci_oc lasttime 0;
   IF FIRST.&outvar THEN DO; 
    IF surv_y<=.z THEN surv_y=lastsurv_y;
    ci_ac=1-surv_y;
  END;
   ci_lc  + MAX((lastsurv_y-surv_y)*(d_multi=1), 0);
   ci_oc  + MAX((lastsurv_y-surv_y)*(d_multi=2), 0);
   *ci_lw  = MAX(0, 1-surv_l, ci_lw);
   altci_ac = SUM(ci_lc, ci_oc);
   IF surv_y<=.z THEN surv_y = 1-ci_ac;
  OUTPUT test; 
  IF LAST.&OUTVAR THEN DO;
     IF  FLOOR(&outvar / 10) NE  FLOOR(lasttime / 10) THEN mark=1;
     ELSE mark = 0;
     OUTPUT gformula.cidata_&outvar._&ds&suff;
     lasttime = &outvar;
     lastsurv_y=surv_y;
  END;
 
 RUN;
 PROC PRINT DATA= gformula.cidata_&outvar._&ds&suff(WHERE=(mark=1));
 RUN;
%MEND;

*make a smaller dataset with start/stop times (rather than person-period data) to reduce computation for phreg;
PROC SQL;
 CREATE TABLE gf AS SELECT 
   intervention, gfid, study, d_lc, d_other, sampwt,
   MIN(agein) AS agest, MAX(age) AS ageend 
  FROM gformula
  GROUP BY gfid
  HAVING age = ageend
  ORDER BY intervention, gfid, ageend;
QUIT;
PROC SQL; DROP TABLE gformula, bootsample, an; QUIT;


*INTERVENTION 0: Natural course;
%LET int=0;
%mksurvdata(DS=gf, invar=agest, outvar=ageend, suff=&int, intervention=&int);
*INTERVENTION 1: dynamic intervention, occ limit < 100 mg/m^3;
*limit = 0.1;
%LET int=1;
%mksurvdata(DS=gf, invar=agest, outvar=ageend, suff=&int, intervention=&int);
*INTERVENTION 2: dynamic intervention, occ limit < 50 mg/m^3;
*limit = 0.05;
%LET int=2;
%mksurvdata(DS=gf, invar=agest, outvar=ageend, suff=&int, intervention=&int);
*INTERVENTION 3: dynamic intervention, occ limit < 25 mg/m^3;
*limit = 0.025;
%LET int=3;
%mksurvdata(DS=gf, invar=agest, outvar=ageend, suff=&int, intervention=&int);

*INTERVENTION 4: dynamic intervention, occ limit < 10 mg/m^3;
*limit = 0.01;
%LET int=4;
%mksurvdata(DS=gf, invar=agest, outvar=ageend, suff=&int, intervention=&int);

**INTERVENTION -1: never exposed;
*limit = 0.00;
%LET int=-1;
%mksurvdata(DS=gf, invar=agest, outvar=ageend, suff=m1, intervention=&int);

RUN;QUIT;RUN;

