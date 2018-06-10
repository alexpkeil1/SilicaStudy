/**********************************************************************************************************************
* Author: Alex Keil
* Project: Pooled silica exposed worker study
* Released under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html
* this program carries out modeling and Monte Carlo algorithms for estimating the parameters of the g-formula
**********************************************************************************************************************/
LIBNAME an  "data/an";


%INCLUDE "transformations.sas"; *variable transformations;
%LET nmc= 300000; * number of Monte Carlo samples;


*step 0 - data preparation;
DATA an;
 SET an.an0001;
 WHERE age<91;
 *SET an.an0001_samp;
 d_other = d_any-d_lc;
 s1 = (study=1);
 s2 = (study=2);
 s3 = (study=3);
 s4 = (study=4);
 s5 = (study=5);
 s6 = (study=6);
 s7 = (study=7);
 s8 = (study=8);
 s9 = (study=9);
 s10 = (study=10);

IF study IN (1,3,5,7) THEN country=1; *us;
 ELSE IF study IN (8,9,10) THEN country=2; *china;
 ELSE country=3; *finns, australians;
 c1 = (country=1);
 c2 = (country=2);

IF age >= ageeof THEN censadmin=1; ELSE censadmin=0;
 IF exprt>0 THEN logexprt=LOG(exprt);
 ELSE logexprt=.z;
 %TRANSFORM();
 %TRANSFORMEXP();
 IF study IN (8,9,10) THEN threshold = -5;
 ELSE IF study IN (1) THEN threshold = -8;
 ELSE threshold = -1e12;
 IF yrsexp>0 THEN anyexp = (exp/yrsexp > EXP(threshold));
 ELSE anyexp = 0;
 IF yrsexplag1>0 THEN anyexplag1 = (explag1/yrsexplag1 > EXP(threshold));
 ELSE anyexplag1 = 0;
RUN;

PROC SORT DATA = an;
 BY study id agein;
 
*create baseline data for MC sampling;
DATA anfirst;
 SET an;
 BY study id agein;
 IF first.id THEN OUTPUT;
RUN;

************************************************************************;
*step 1 - modeling;
************************************************************************;


*exposure, any;
%LET aapreds = agec agesq agecu datec datesq datecu
               sex race
               exp_1_15
               cye1_3_1 cye1_3_2;

*exposure, how much (Non-chinese);
%LET apreds =  agec agework_sp1 datec datework_sp1 datework_sp2
               sex race
               explag1 explag1*explag1 exp_2_15 exp_2_15*exp_2_15 cumexplag15 cumexplag15*cumexplag15
               cumyrsexplag1 cyxl1sq cyxl1cu;

*exposure, how much (Chinese only);
%LET apredsC = agec ageworkchina_sp1 datec dateworkchina_sp1 dateworkchina_sp2
               sex
               xl1cat_5_2 xl1cat_5_3 xl1cat_5_4 sil2_4_1 sil2_4_2 sil2_4_3
               cumyrsexplag1 cyxl1sq cyxl1cu;

*leaving employment;
%LET lpreds = agec agesq agecu datec datesq datecu
              sex race
              cumexplag1;

*lung cancer;
%LET y1preds = agec agedeath3b_sp1 agedeath3b_sp2 agedeath3b_sp3
               datec datesq datecu
               sex /*race*/
               atworklag1 
               exp_1_15 exp_1_15*exp_1_15  cumexplag15 cumexplag15*cumexplag15 
               cumyrsexplag1 cyxl1sq cyxl1cu
               ;

*lung cancer - Chinese cohorts - need different baseline risk function;
%LET y1predsC = agec agedeathchina4b_sp1 agedeathchina4b_sp2 agedeathchina4b_sp3 agedeathchina4b_sp4
                datec datesq datecu
                sex
                atworklag1 
                exp_1_15 exp_1_15*exp_1_15  cumexplag15 cumexplag15*cumexplag15 
                cumyrsexplag1 cyxl1sq cyxl1cu
                ;


*other causes;
%LET y2preds = agec agedeath3a_sp1 agedeath3a_sp2 agedeath3a_sp3
               datec datesq datecu
               sex /*race*/
               atworklag1
               exp_1_15 exp_1_15*exp_1_15  cumexplag15 cumexplag15*cumexplag15 
               cumyrsexplag1 cyxl1sq cyxl1cu
               ;

*other causes - Chinese cohorts;
%LET y2predsC = agec agedeathchina4a_sp1 agedeathchina4a_sp2 agedeathchina4a_sp3 agedeathchina4a_sp4
                datec datesq datecu
                sex
                atworklag1 
                exp_1_15 exp_1_15*exp_1_15  cumexplag15 cumexplag15*cumexplag15 
                cumyrsexplag1 cyxl1sq cyxl1cu
                ;

*censoring;
%LET cpreds =  agec agedeath3a_sp1 agedeath3a_sp2 agedeath3a_sp3
               datec datesq datecu
               sex race
               atworklag1
               exp_1_15 exp_1_15*exp_1_15  cumexplag15 cumexplag15*cumexplag15 
               cumyrsexplag1 cyxl1sq cyxl1cu
               ;

%LET censstudies = 1,3,4,7; *studies with censoring;
%LET cstudies = 8,9,10; *chinese cohorts - very different baseline risk;
%LET noncstudies = 1,2,3,4,5,6,7; * non chinese cohorts;
%LET zerostudies = 1,8,9,10; *studies with some unexposed workers;

/*************************************** modeling ************************************************/

*exposure - any;
PROC GENMOD DATA = an DESCENDING ;
 BY study;
 WHERE study IN (&zerostudies) AND yrsexp>0; *whether any exposure occurred while at work!;
 MODEL anyexp = &aaPREDS
  / D=B LINK=LOGIT;
 ODS OUTPUT ParameterEstimates=Aa_model (KEEP=study parameter estimate  WHERE=(parameter ^= "Scale"));
RUN;

*exposure - quantitative | any;
PROC GENMOD DATA = an ;
BY study;
 WHERE STUDY IN(&noncstudies) AND yrsexp>0 AND anyexp=1; *modeling exposure rate (will be one year per record in simulated data);
 MODEL logexprt = &aPREDS
  / D=NORMAL LINK=IDENTITY;
 ODS OUTPUT ParameterEstimates=A_model (KEEP=study parameter estimate );
RUN;

*Chinese cohorts;
PROC GENMOD DATA = an ;
BY study;
 WHERE STUDY IN (&cstudies) AND yrsexp>0 AND anyexp=1 ; *modeling exposure rate (will be one year per record in simulated data);
 MODEL logexprt = &aPREDSC
  / D=NORMAL LINK=IDENTITY;
 ODS OUTPUT ParameterEstimates=A_modelC (KEEP=study parameter estimate );
RUN;

*employment status;
PROC GENMOD DATA = an DESCENDING;
BY study;
 WHERE yrsexp>0 OR leftwork=1 ;
 MODEL leftwork = &lPREDS
  / D=B LINK=LOGIT;
 ODS OUTPUT ParameterEstimates=L_model (KEEP=study parameter estimate  WHERE=(parameter ^= "Scale"));
RUN;

*lung cancer deaths ;
PROC GENMOD DATA = an DESCENDING;
WHERE STUDY IN (&noncstudies) ;
BY study;
 MODEL d_lc = &y1preds
  / D=B LINK=LOGIT;
 ODS OUTPUT ParameterEstimates=Y1_model (KEEP=study parameter estimate  WHERE=(parameter ^= "Scale"));
RUN;

*lung cancer deaths - Chinese cohorts ;
PROC GENMOD DATA = an DESCENDING;
 WHERE study IN (&cstudies);
 BY study;
 MODEL d_lc = &y1predsc
  / D=B LINK=LOGIT;
 ODS OUTPUT ParameterEstimates=Y1_modelc (KEEP=study parameter estimate  WHERE=(parameter ^= "Scale"));
RUN;

*deaths attributed to other causes ;
PROC GENMOD DATA = an DESCENDING;
WHERE STUDY IN (&noncstudies) ;
BY study;
 MODEL d_other = &y2preds
  / D=B LINK=LOGIT;
 ODS OUTPUT ParameterEstimates=Y2_model (KEEP=study parameter estimate WHERE=(parameter ^= "Scale"));
RUN;

*deaths attributed to other causes - Chinese cohorts ;
PROC GENMOD DATA = an DESCENDING;
WHERE study IN (&cstudies) ;
BY study;
 MODEL d_other = &y2predsc
  / D=B LINK=LOGIT;
 ODS OUTPUT ParameterEstimates=Y2_modelc (KEEP=study parameter estimate WHERE=(parameter ^= "Scale"));
RUN;

*censoring;
PROC GENMOD DATA = an DESCENDING;
WHERE STUDY IN (&censstudies) ;
BY study;
 MODEL ltfu = &cpreds
  / D=B LINK=LOGIT;
 ODS OUTPUT ParameterEstimates=C_model (KEEP=study parameter estimate  WHERE=(parameter ^= "Scale"));
RUN;


/************************************************************************************************************************
****** end modeling
************************************************************************************************************************/

*max exposure;
PROC MEANS DATA = an NOPRINT;
 BY study;
 VAR exprt;
 OUTPUT OUT=maxexp(DROP=_:) MAX=max_exp MIN=min_exp;
RUN;

************************************************************************;
*step 2 - MC sample of baseline, macro automation of covariates;
************************************************************************;
PROC FREQ DATA = an NOPRINT; TABLES study / OUT=studyn(KEEP=study percent); RUN;
DATA studyn; SET studyn; _NSIZE_=ROUND(.01*percent*&nmc);RUN;

*equal allocation to each dataset (requires weights to get back sample);
PROC SURVEYSELECT DATA = anfirst OUT=mcsample SAMPSIZE=&NMC METHOD=URS OUTHITS;
 STRATA study / ALLOC=(0.0796 0.0796 0.0796 0.0796 0.0796 0.0796 0.0796 0.1476 0.1476 0.1476); * oversample Chinese cohorts due to low lung cancer risk;
RUN;

PROC FREQ DATA = anfirst NOPRINT; TABLES study / OUT = anwts; RUN;
PROC FREQ DATA = mcsample NOPRINT; TABLES study / OUT = anwts2; RUN; *inefficient, but automatic way of just getting the allocation sample proportions back from surveyselect;
*weights to re-weight sample;
DATA sampwt(KEEP=study sampwt);
 SET anwts;
 sampwt = percent/100; 
RUN;
DATA w1(KEEP=study w1);
 SET anwts;
 w1 = percent/100; 
RUN;
DATA w2(KEEP=study w2);
 SET anwts2;
 w2 = percent/100; 
RUN;

DATA sampwt(KEEP=study sampwt w1 w2);
 MERGE w1 w2;
 sampwt=w1/w2;
RUN;
PROC PRINT DATA = sampwt; TITLE "Allocation by study in cohort (w1) into MC sample (w2) and analytic weights (sampwt)";RUN;


%MACRO count();
* counting coefficients from each model;
%GLOBAL aamod  amod amodc lmod y1mod y2mod y1modc y2modc cmod
        naa na nac namc1 namc2  nmmc1 nmmc2  nl ny1 ny2 ny1c ny2c nc;
/* exposure */
%LET i=1; %LET amod=_a1; %DO %UNTIL(%QSCAN(&apreds, %EVAL(&i), " ")=);
   %LET amod=&amod + _a%EVAL(&i+1) * %QSCAN(&apreds, %EVAL(&i), " ");  %LET i = %EVAL(&i+1);
%END;%LET na = &i;
%LET i=1; %LET amodc=_ac1; %DO %UNTIL(%QSCAN(&apredsc, %EVAL(&i), " ")=);
   %LET amodc=&amodc + _ac%EVAL(&i+1) * %QSCAN(&apredsc, %EVAL(&i), " ");  %LET i = %EVAL(&i+1);
%END;%LET nac = &i;

%LET i=1; %LET aamod=_aa1; %DO %UNTIL(%QSCAN(&aapreds, %EVAL(&i), " ")=);
   %LET aamod=&aamod + _aa%EVAL(&i+1) * %QSCAN(&aapreds, %EVAL(&i), " ");  %LET i = %EVAL(&i+1);
%END;%LET naa = &i;
/* employment status */
%LET i=1; %LET lmod=_l1; ;%DO %UNTIL(%QSCAN(&lpreds, %EVAL(&i), " ")=);
   %LET lmod=&lmod + _l%EVAL(&i+1) * %QSCAN(&lpreds, %EVAL(&i), " "); %LET i = %EVAL(&i+1);
 %END;%LET nl = &i;

/* lung cancer outcomes */
%LET i=1; %LET y1mod=_da1; %DO %UNTIL(%QSCAN(&y1preds, %EVAL(&i), " ")=);
   %LET y1mod=&y1mod + _da%EVAL(&i+1) * %QSCAN(&y1preds, %EVAL(&i), " ");  %LET i = %EVAL(&i+1);
%END;%LET ny1 = &i;
%LET i=1; %LET y1modc=_dac1; %DO %UNTIL(%QSCAN(&y1predsc, %EVAL(&i), " ")=);
   %LET y1modc=&y1modc + _dac%EVAL(&i+1) * %QSCAN(&y1predsc, %EVAL(&i), " ");  %LET i = %EVAL(&i+1);
%END;%LET ny1c = &i;

/* other mortality outcomes */
%LET i=1; %LET y2mod=_db1; %DO %UNTIL(%QSCAN(&y2preds, %EVAL(&i), " ")=);
   %LET y2mod=&y2mod + _db%EVAL(&i+1) * %QSCAN(&y2preds, %EVAL(&i), " ");  %LET i = %EVAL(&i+1);
%END;%LET ny2 = &i;
%LET i=1; %LET y2modc=_dbc1; %DO %UNTIL(%QSCAN(&y2predsc, %EVAL(&i), " ")=);
   %LET y2modc=&y2modc + _dbc%EVAL(&i+1) * %QSCAN(&y2predsc, %EVAL(&i), " ");  %LET i = %EVAL(&i+1);
%END;%LET ny2c = &i;
/* censoring */
%LET i=1; %LET cmod=_c1; %DO %UNTIL(%QSCAN(&cpreds, %EVAL(&i), " ")=);
   %LET cmod=&cmod + _c%EVAL(&i+1) * %QSCAN(&cpreds, %EVAL(&i), " ");  %LET i = %EVAL(&i+1);
%END;%LET nc = &i;
%MEND;
%COUNT;
*check log for what the data generating models look like;

* process coefficients for use in mc sampling;
*mortality;
PROC TRANSPOSE DATA = Y1_model OUT=c_d1(DROP=_NAME_); BY STUDY; ID parameter; ;RUN;
PROC TRANSPOSE DATA = Y2_model OUT=c_d2(DROP=_NAME_); BY STUDY; ID parameter; ;RUN;
PROC TRANSPOSE DATA = Y1_modelC OUT=c_d1c(DROP=_NAME_); BY STUDY; ID parameter; ;RUN;
PROC TRANSPOSE DATA = Y2_modelC OUT=c_d2c(DROP=_NAME_); BY STUDY; ID parameter; ;RUN;
*employment;
PROC TRANSPOSE DATA = L_model OUT=c_l(DROP=_NAME_); BY STUDY; ID parameter; ;RUN;
*exposure;
PROC TRANSPOSE DATA = A_model OUT=c_a(DROP=_NAME_); BY STUDY; ID parameter; ;RUN;
PROC TRANSPOSE DATA = Aa_model OUT=c_aa(DROP=_NAME_); BY STUDY; ID parameter; ;RUN;
PROC TRANSPOSE DATA = A_modelc OUT=c_ac(DROP=_NAME_); BY STUDY; ID parameter; ;RUN;
*censoring;
PROC TRANSPOSE DATA = C_model OUT=c_c(DROP=_NAME_); BY STUDY; ID parameter; ;RUN;

DATA st(KEEP=study); SET c_l;
DATA sta(KEEP=study); SET c_aa;
DATA stc(KEEP=study); SET c_ac;RUN;
DATA stcens(KEEP=study); SET c_c;RUN;

DATA c_l; SET c_l(DROP=study);
 ARRAY coefs[*] _NUMERIC_; ARRAY _l[&nl];
 DO j = 1 TO DIM(coefs); _l[j] = coefs[j];END;
 KEEP _l:;
DATA c_d1; SET c_d1(DROP=study);
 ARRAY coefs[*] _NUMERIC_;ARRAY _da[&ny1];
 DO j = 1 TO DIM(coefs);_da[j] = coefs[j];END;
 KEEP _da:;
DATA c_d2; SET c_d2(DROP=study);
 ARRAY coefs[*] _NUMERIC_;ARRAY _db[&ny2];
 DO j = 1 TO DIM(coefs);_db[j] = coefs[j];END;
 KEEP _db:;
DATA c_a; SET c_a(DROP=study);
 ARRAY coefs[*] _NUMERIC_;ARRAY _a[&na];
 DO j = 1 TO (DIM(coefs)-1); _a[j] = coefs[j];END;
 _a_scale = scale;
 KEEP _a:;
DATA c_aa; SET c_aa(DROP=study);
 ARRAY coefs[*] _NUMERIC_;ARRAY _aa[&naa];
 DO j = 1 TO DIM(coefs); _aa[j] = coefs[j];END;
 KEEP _aa:;
RUN;
DATA c_ac; SET c_ac(DROP=study);
 ARRAY coefs[*] _NUMERIC_;ARRAY _ac[&nac];
 DO j = 1 TO (DIM(coefs)-1); _ac[j] = coefs[j];END;
 _ac_scale = scale;
 KEEP _ac:;
DATA c_d1c; SET c_d1c(DROP=study);
 ARRAY coefs[*] _NUMERIC_;ARRAY _dac[&ny1c];
 DO j = 1 TO DIM(coefs);_dac[j] = coefs[j];END;
 KEEP _dac:;
DATA c_d2c; SET c_d2c(DROP=study);
 ARRAY coefs[*] _NUMERIC_;ARRAY _dbc[&ny2c];
 DO j = 1 TO DIM(coefs);_dbc[j] = coefs[j];END;
 KEEP _dbc:;
DATA c_c; SET c_c(DROP=study);
 ARRAY coefs[*] _NUMERIC_;ARRAY _c[&nc];
 DO j = 1 TO DIM(coefs);_c[j] = coefs[j];END;
 KEEP _c:;RUN;
OPTIONS MERGENOBY=NOWARN;
DATA c_l; MERGE st c_l;
DATA c_d1; MERGE st c_d1;
DATA c_d2; MERGE st c_d2;
DATA c_a; MERGE st c_a;
DATA c_aa; MERGE sta c_aa;
DATA c_ac; MERGE stc c_ac;
DATA c_c; MERGE stcens c_c;
*mortality;
DATA c_d1c; MERGE stc c_d1c;
DATA c_d2c; MERGE stc c_d2c;
DATA c_aa; MERGE st c_aa; BY study; RUN;
OPTIONS MERGENOBY=WARN;


************************************************************************;
*step 3 - MC sample of time-varying quantities;
************************************************************************;
DATA gformula;
 MERGE mcsample (KEEP = study id agein datein race sex decob ageeof dob yrsexplag: explag: cumexplag1 cumexplag2 cumexplag15 cumyrsexplag1 cumyrsexplag15)
  c_a c_aa c_ac  c_l  c_d1  c_d2 c_d1c  c_d2c c_c maxexp sampwt;
  BY study;
 newid = _N_;
RUN;
DATA gformula(KEEP=study sampwt gfid intervention agein age datein date race sex decob ageeof dob exp logexp cumexp atwork anyexp yrsexp leftwork cumexplag15 cumyrsexplag15 d_: ltfu censadmin py
              %SYSFUNC(TRANSTRN( &Lpreds &apreds &y1preds &y2preds, * ,%STR( ))));
 LENGTH study gfid newid intervention agein age d_lc d_other d_any ltfu censadmin leftwork exp cumexp atwork anyexp cumyrsexp 8;
 SET gformula;
 FORMAT intervention intform.;

 ARRAY llag[15] yrsexplag1-yrsexplag15;
 ARRAY alag[15] explag1-explag15;
 ARRAY ollag[15] oldyrsexplag1-oldyrsexplag15;
 ARRAY oalag[15] oldexplag1-oldexplag15;
 FORMAT date datein MMDDYY10.;
 oldagein = agein; *keep this for each ID;
 olddatein = datein;
 *exposure;
 oldcumexp  =        cumexp;
 oldcumyrsexp =      cumyrsexp;
 
  DO i = 1 TO 15;
   ollag[i] = llag[i];
   oalag[i] = alag[i];
  END;

 
 oldcumexplag1 =     cumexplag1 ;
 oldcumexplag2 =     cumexplag2 ;
 oldcumexplag15 =    cumexplag15 ;
 oldcumyrsexplag1 =  cumyrsexplag1;
 oldcumyrsexplag15 = cumyrsexplag15;
 py=1;

 *MAIN intervention loop - repeat for each new occupational standard;
 DO intervention = -1, 0, 1, 2, 3, 4;
  gfid = intervention*&nmc+newid+1;

  *initializing time-varying variables;
  agein = (oldagein);
  datein = olddatein;
  *exposure;
  exp = 0;
  *employment;
  atwork=1;
  leftwork=0;
  *years of exposure;
  yrsexp = 1;

  *lagged;
  DO i = 1 TO 15;
   llag[i] = ollag[i];
   alag[i] = oalag[i];
  END;

  cumexplag1 =     oldcumexplag1 ;
  cumexplag2 =     oldcumexplag2 ;
  cumexplag15 =    oldcumexplag15 ;
  cumyrsexplag1 =  oldcumyrsexplag1;
  cumyrsexplag15 = oldcumyrsexplag15;

  firstobs = 1;
  done = 0;
  *mortality and follow up;
  ltfu = 0;
  d_lc = 0;
  d_other = 0;
  *MAIN time loop;
  DO WHILE(done=0);
   *external time varying quantities;
   age=agein+1;
   IF MONTH(datein)=2 AND DAY(datein)=29 THEN DO; *leap years;
    date = MDY(MONTH(datein), 28, YEAR(datein)+1)-1;
   END;
   ELSE DO; *all other dates;
    date = MDY(MONTH(datein), DAY(datein), YEAR(datein)+1)-1;
   END;
   * all variable transformations here;
   %TRANSFORM(); * transform age, date;
   %TRANSFORMEXP(); *transform exposure variables (all lagged);

   ****** internal (modeled) time varying quantities;
    **************;
    /*  employment */
    **************;
   IF firstobs THEN DO;
     firstobs = 0; *do nothing here - every individual employed in first observation;
     yrsexp = 1;
   END;
   ELSE DO;
    IF llag[1]>0 THEN DO;
     *if already at work, simulate probability of leaving work;
     lol = &LMOD;
     IF .z < lol < -700 THEN lol=-700;
     IF lol > 700  THEN lol=700;
     atwork = 1 - RAND('BERNOULLI', 1/(1+exp(-lol)));
     yrsexp = atwork;
     IF atwork=0 THEN leftwork=1;
     ELSE leftwork = 0;
    END; * llag[1]>0;
    ELSE DO;
     *if not already at work, do nothing;
     atwork = 0;
     yrsexp = 0;
     leftwork = 0;
    END;
   END;
    **************;
    /*  exposure */
    **************;
   IF atwork>0 THEN DO;
    *some studies have unexposed - this is the log odds of any exposure (hurdle model);
    IF study IN (&zerostudies) THEN DO;
     loa = &aaMOD;
     IF .z < loa < -700 THEN loa=-700;
     IF loa > 700  THEN loa=700;
    END; *study IN zerostudies;

    *studies with assumed normal model for quantitative exposure;
    IF study IN (&noncstudies) THEN DO;
      ma = &AMOD;
      _SIGMA = _a_scale;*standard deviation;
    END; *study IN noncstudies;

    *modeling exposure among the c-studies;
    IF study IN (&cstudies) THEN DO;
      ma = &amodc;
      _SIGMA = _ac_SCALE; *standard deviation;
    END; *study in cstudies;

    /* BEGIN INTERVENTIONS */
    /* natural course or interventions involving the natural course */
    IF intervention >= 0 OR intervention = -2 THEN DO;
     *some studies had unexposed individuals - simulate that here;
     IF study IN (&zerostudies) THEN anyexp = RAND('BERNOULLI', 1/(1+exp(-loa)));
     ELSE anyexp=1;
     *simulate quantitative exposure among those exposed here - truncated at observed max/min;
     IF anyexp = 1 THEN DO;
      count=1;
      countmin=1;
      ******* main exposure generating statement *********;
      makeexp: logexp = RAND('NORMAL', ma, _SIGMA); *mean, sd;
      *truncating extreme exposures - loop 40 times resampling from the distribution, otherwise just take the limit;
      /* natural course */
      IF intervention  IN (0, -2) THEN DO;
       IF EXP(logexp) > max_exp OR EXP(logexp) < min_exp THEN DO;
        count=count+1;
        IF count<=40 THEN GOTO makeexp;
        ELSE IF count>40 THEN DO;
         IF logexp > LOG(max_exp) THEN logexp=LOG(max_exp);
         ELSE IF logexp < LOG(min_exp) THEN logexp=LOG(min_exp);
        END; *count>40;
       END; *EXP(logexp) > max_exp;
      END; *intervention 0;

      /* occupational limit 1 */
      IF intervention  = 1 THEN DO;
       limit = 0.1;
       *intervention: old PEL = 100 ug/m^3 = 0.1;
       IF logexp > LOG(limit) THEN logexp = LOG(limit);
       *keep lower bound in the observed range;
       IF EXP(logexp) < min_exp THEN DO;
        count=count+1;
        IF count<=40 THEN GOTO makeexp;
        ELSE logexp=LOG(min_exp);
       END; *EXP(logexp) < min_exp;
      END; *intervention 1;

      /* occupational limit 2 */
      IF intervention  = 2 THEN DO;
       limit = 0.05;
       *intervention: new PEL = 50 ug/m^3 = 0.05;
       IF logexp > LOG(limit) THEN logexp = LOG(limit);
       *keep lower bound in the observed range;
       IF EXP(logexp) < min_exp THEN DO;
        count=count+1;
        IF count<=40 THEN GOTO makeexp;
        ELSE logexp=LOG(min_exp);
       END; *EXP(logexp) < min_exp;
      END; *intervention 2;

      /* occupational limit 3 */
      IF intervention  = 3 THEN DO;
       limit = 0.025;
       *intervention: alt PEL = 25 ug/m^3 = 0.01;
       IF logexp > LOG(limit) THEN logexp = LOG(limit);
       *keep lower bound in the observed range;
       IF EXP(logexp) < min_exp THEN DO;
        count=count+1;
        IF count<=40 THEN GOTO makeexp;
        ELSE logexp=LOG(min_exp);
       END; *EXP(logexp) < min_exp;
      END; *intervention 3;

      /* occupational limit 4 */
      IF intervention  = 4 THEN DO;
       limit = 0.01;
       *intervention: alt PEL = 10 ug/m^3 = 0.01;
       IF logexp > LOG(limit) THEN logexp = LOG(limit);
       *keep lower bound in the observed range;
       IF EXP(logexp) < min_exp THEN DO;
        count=count+1;
        IF count<=40 THEN GOTO makeexp;
        ELSE logexp=LOG(min_exp);
       END; *EXP(logexp) < min_exp;
      END; *intervention 3;


      /* END RANDOM DYNAMIC INTERVENTIONS BASED ON NATURAL COURSE */
     END; *anyexp = 1;
     ELSE IF anyexp = 0 THEN DO;
      *the level of log exposure rate assumed among those unexposed;
      IF study IN (&cstudies) THEN logexp = -7;
      ELSE IF study IN (1) THEN logexp = -11;
     END; *anyexp = 0;
     *exponentiated exposure;
     exp = EXP(logexp);
    END; *interventions based on natural course (>=0);

    /* DETERMINISTIC INTERVENTIONS */
    ELSE IF intervention=-1 THEN DO;
     * never exposed;
	 anyexp=0;
     exp = 0;
     logexp = LOG(1/365);
    END; *intervention=-1;
    /* END DETERMINISTIC INTERVENTIONS */
   END;*atwork;
    **************;
    /*  end employment */
    **************;
   IF atwork=0 THEN DO;
    anyexp=0;
    exp=0;
    yrsexp=0;
    logexp = LOG(1/365);
   END;

   *cumulative variables;
   cumexp = SUM(cumexplag1, exp);
   cumyrsexp = SUM(cumyrsexplag1, yrsexp);

   *mortality - the predictors vary by study;
   IF study IN (&noncstudies) THEN lod1 = &y1mod;
   IF study IN (&cstudies) THEN lod1 = &y1modc;
   IF lod1 < -700 THEN lod1 = -700;
   IF lod1 > 700  THEN lod = 700;
   d_lc = RAND('BERNOULLI', 1/(1+exp(-lod1)));

   IF d_lc = 0 THEN DO;
    IF study IN(&noncstudies) THEN lod2 = &y2mod;
    IF study IN(&cstudies) THEN lod2 = &y2modc;
    IF lod2 < -700 THEN lod2 = -700;
    IF lod2 > 700  THEN lod2 = 700;
    d_other = RAND('BERNOULLI', 1/(1+exp(-lod2)));
   END; *d_lc = 0;
   ELSE d_other=0;

   d_any = d_lc+d_other;
   *censoring by loss to follow-up;
   IF intervention=-2 THEN DO; *natural course except censoring;
    IF d_any = 0 AND study IN (&censstudies) THEN DO;
     loc = &cmod;
     IF .z < loc < -700 THEN loc=-700;
     IF loc > 700  THEN loc=700;
     ltfu = RAND('BERNOULLI', 1/(1+exp(-loc)));
    END;
   END;

   *exit if reach age 90 or censoring date;
   IF age >= ageeof AND d_any=0 THEN censadmin=1; ELSE censadmin=0;
   IF age >= 90 OR d_any=1 OR censadmin=1 OR ltfu=1 THEN done=1;
   *dont include those younger than 16 years old;
   IF agein >= 16 THEN OUTPUT;

   *lagged exposure and employment variables;
   DO i = 15 TO 2 BY -1;
    llag[i] = llag[i-1];
    alag[i] = alag[i-1];
   END;
   llag[1] = yrsexp;
   alag[1] = exp;
   cumexplag1 = SUM(cumexplag1, alag[1]);
   cumexplag2 = SUM(cumexplag2, alag[2]);
   cumexplag5 = SUM(cumexplag5, alag[5]);
   cumexplag15 = SUM(cumexplag15, alag[15]);
   cumyrsexplag1 = SUM(cumyrsexplag1, llag[1]);
   cumyrsexplag15 = SUM(cumyrsexplag15, llag[15]);
   IF cumexplag15>0 THEN logcumexplag15 = LOG(cumexplag15); ELSE LOGcumexplag15 = LOG(1/365);
   *incrementing age;
   agein=agein+1; *increment age by one year from last entry;
   datein=date+1; *increment date by one day from last exit;
   anyexplag1=anyexp;
   yrsexplag1=yrsexp;
   explag1=exp;
  END;*time loop;
 END; *intervention loop;
RUN;

