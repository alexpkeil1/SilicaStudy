%MACRO transform();
    *transformations for both the observed and g-formula data, external time variables;
    *age;
    agec=agein-45;
    agesq=agec*agec;
    agecu=agesq*agec;

    *date;
    datec=datein/5000;
    datesq=datec*datec;
    datecu=datesq*datec;

    *splines for age - based on quantiles given by Frank Harrell;
    * different knots for at work vs. all person time, different causes of death, also different for Chinese cohorts due to very different age/mortality patterns;
    *knots: -22 -3 18 ages: 22, 42, 63;
    agework_sp1 = ((agec>-22) * (agec--22)/((18--22)**(2/3)))**3
        + ((-3--22)* ((agec>18) * (agec-18)/((18--22)**(2/3)))**3 - 
           (18--22) * ((agec>-3) * (agec--3)/((18--22)**(2/3)))**3)/(18--3);
        
    *knots: -0.74 0.168 0.821 1.842 1949-11-14 1962-04-20 1971-03-29 1985-03-20;
    datework_sp1 = ((datec>-0.74) * (datec--0.74)/((1.842--0.74)**(2/3)))**3
        + ((0.821--0.74)* ((datec>1.842) * (datec-1.842)/((1.842--0.74)**(2/3)))**3 - 
           (1.842--0.74) * ((datec>0.821) * (datec-0.821)/((1.842--0.74)**(2/3)))**3)/(1.842-0.821);

    datework_sp2 = ((datec>0.168) * (datec-0.168)/((1.842--0.74)**(2/3)))**3
        + ((0.821-0.168)* ((datec>1.842) * (datec-1.842)/((1.842--0.74)**(2/3)))**3 - 
           (1.842-0.168) * ((datec>0.821) * (datec-0.821)/((1.842--0.74)**(2/3)))**3)/(1.842-0.821);

    *knots: -20 -3 11 ages: 24, 42, 56;
    ageworkchina_sp1 = ((agec>-20) * (agec--20)/((11--20)**(2/3)))**3
        + ((-3--20)* ((agec>11) * (agec-11)/((11--20)**(2/3)))**3 - 
           (11--20) * ((agec>-3) * (agec--3)/((11--20)**(2/3)))**3)/(11--3);


    *knots: 0.94 1.306 1.744 2.45 1972-11-13 1977-11-17 1983-11-16 1993-07-16;
    dateworkchina_sp1 = ((datec>0.94) * (datec-0.94)/((2.45-0.94)**(2/3)))**3
        + ((1.744-0.94)* ((datec>2.45) * (datec-2.45)/((2.45-0.94)**(2/3)))**3 - 
           (2.45-0.94) * ((datec>1.744) * (datec-1.744)/((2.45-0.94)**(2/3)))**3)/(2.45-1.744);

    dateworkchina_sp2 = ((datec>1.306) * (datec-1.306)/((2.45-0.94)**(2/3)))**3
        + ((1.744-1.306)* ((datec>2.45) * (datec-2.45)/((2.45-0.94)**(2/3)))**3 - 
           (2.45-1.306) * ((datec>1.744) * (datec-1.744)/((2.45-0.94)**(2/3)))**3)/(2.45-1.744);

    *knots: -1.4 12.1 18.2 29.1 - 43.6 57.1 63.2 74.1;
    *agedeathchina;
    *knots: 1.136 1.872 2.206 2.462 [1] 1975-07-20 1985-08-17 1990-03-14 1993-09-14;
    datedeathchina_sp1 = ((datec>1.136)*((datec-1.136)/1.136)**3) 
        +((datec>2.462)*((datec-2.462)/1.136)**3)*(2.206-1.136) 
        -((datec>2.206)*((datec-2.206)/1.136)**3)*(2.462-1.136)/(2.462-2.206);

    datedeathchina_sp2 = ((datec>1.872)*((datec-1.872)/1.136)**3) 
        +((datec>2.462)*((datec-2.462)/1.136)**3)*(2.206-1.872) 
        -((datec>2.206)*((datec-2.206)/1.136)**3)*(2.462-1.872)/(2.462-2.206);
        
    * different knots for other causes due to later times of death than lung cancer;
    *knots: -3.1 13 19.8 26.9 38.3  - non lung cancer;
    agedeath3a_sp1 = ((agec>-3.1) * (agec--3.1)/((38.3--3.1)**(2/3)))**3
        + ((26.9--3.1)* ((agec>38.3) * (agec-38.3)/((38.3--3.1)**(2/3)))**3 - 
           (38.3--3.1) * ((agec>26.9) * (agec-26.9)/((38.3--3.1)**(2/3)))**3)/(38.3-26.9);

    agedeath3a_sp2 = ((agec>13) * (agec-13)/((38.3--3.1)**(2/3)))**3
        + ((26.9-13)* ((agec>38.3) * (agec-38.3)/((38.3--3.1)**(2/3)))**3 - 
           (38.3-13) * ((agec>26.9) * (agec-26.9)/((38.3--3.1)**(2/3)))**3)/(38.3-26.9);

    agedeath3a_sp3 = ((agec>19.8) * (agec-19.8)/((38.3--3.1)**(2/3)))**3
        + ((26.9-19.8)* ((agec>38.3) * (agec-38.3)/((38.3--3.1)**(2/3)))**3 - 
           (38.3-19.8) * ((agec>26.9) * (agec-26.9)/((38.3--3.1)**(2/3)))**3)/(38.3-26.9);

    *knots: 4.8 15 20.8 26 35.1 - lung cancer;
    agedeath3b_sp1 = ((agec>4.8) * (agec-4.8)/((35.1-4.8)**(2/3)))**3
        + ((26-4.8)* ((agec>35.1) * (agec-35.1)/((35.1-4.8)**(2/3)))**3 - 
           (35.1-4.8) * ((agec>26) * (agec-26)/((35.1-4.8)**(2/3)))**3)/(35.1-26);

    agedeath3b_sp2 = ((agec>15) * (agec-15)/((35.1-4.8)**(2/3)))**3
        + ((26-15)* ((agec>35.1) * (agec-35.1)/((35.1-4.8)**(2/3)))**3 - 
           (35.1-15) * ((agec>26) * (agec-26)/((35.1-4.8)**(2/3)))**3)/(35.1-26);

    agedeath3b_sp3 = ((agec>20.8) * (agec-20.8)/((35.1-4.8)**(2/3)))**3
        + ((26-20.8)* ((agec>35.1) * (agec-35.1)/((35.1-4.8)**(2/3)))**3 - 
           (35.1-20.8) * ((agec>26) * (agec-26)/((35.1-4.8)**(2/3)))**3)/(35.1-26);


    *knots: -7 5.3 11.7 16.9 22.8 31.9 ;
    agedeathchina4a_sp1 = ((agec>-7) * (agec--7)/((31.9--7)**(2/3)))**3
        + ((22.8--7)* ((agec>31.9) * (agec-31.9)/((31.9--7)**(2/3)))**3 - 
           (31.9--7) * ((agec>22.8) * (agec-22.8)/((31.9--7)**(2/3)))**3)/(31.9-22.8);

    agedeathchina4a_sp2 = ((agec>5.3) * (agec-5.3)/((31.9--7)**(2/3)))**3
        + ((22.8-5.3)* ((agec>31.9) * (agec-31.9)/((31.9--7)**(2/3)))**3 - 
           (31.9-5.3) * ((agec>22.8) * (agec-22.8)/((31.9--7)**(2/3)))**3)/(31.9-22.8);

    agedeathchina4a_sp3 = ((agec>11.7) * (agec-11.7)/((31.9--7)**(2/3)))**3
        + ((22.8-11.7)* ((agec>31.9) * (agec-31.9)/((31.9--7)**(2/3)))**3 - 
           (31.9-11.7) * ((agec>22.8) * (agec-22.8)/((31.9--7)**(2/3)))**3)/(31.9-22.8);

    agedeathchina4a_sp4 = ((agec>16.9) * (agec-16.9)/((31.9--7)**(2/3)))**3
        + ((22.8-16.9)* ((agec>31.9) * (agec-31.9)/((31.9--7)**(2/3)))**3 - 
           (31.9-16.9) * ((agec>22.8) * (agec-22.8)/((31.9--7)**(2/3)))**3)/(31.9-22.8);

    *knots: -1.4 9.5 13.4 16.7 21.6 29.1 ;
    agedeathchina4b_sp1 = ((agec>-1.4) * (agec--1.4)/((29.1--1.4)**(2/3)))**3
        + ((21.6--1.4)* ((agec>29.1) * (agec-29.1)/((29.1--1.4)**(2/3)))**3 - 
           (29.1--1.4) * ((agec>21.6) * (agec-21.6)/((29.1--1.4)**(2/3)))**3)/(29.1-21.6);

    agedeathchina4b_sp2 = ((agec>9.5) * (agec-9.5)/((29.1--1.4)**(2/3)))**3
        + ((21.6-9.5)* ((agec>29.1) * (agec-29.1)/((29.1--1.4)**(2/3)))**3 - 
           (29.1-9.5) * ((agec>21.6) * (agec-21.6)/((29.1--1.4)**(2/3)))**3)/(29.1-21.6);

    agedeathchina4b_sp3 = ((agec>13.4) * (agec-13.4)/((29.1--1.4)**(2/3)))**3
        + ((21.6-13.4)* ((agec>29.1) * (agec-29.1)/((29.1--1.4)**(2/3)))**3 - 
           (29.1-13.4) * ((agec>21.6) * (agec-21.6)/((29.1--1.4)**(2/3)))**3)/(29.1-21.6);

    agedeathchina4b_sp4 = ((agec>16.7) * (agec-16.7)/((29.1--1.4)**(2/3)))**3
        + ((21.6-16.7)* ((agec>29.1) * (agec-29.1)/((29.1--1.4)**(2/3)))**3 - 
           (29.1-16.7) * ((agec>21.6) * (agec-21.6)/((29.1--1.4)**(2/3)))**3)/(29.1-21.6); 
%MEND;
%MACRO transformexp();
    *exposure categories in Chinese cohorts;
    xl1cat_5_1 = (0 < explag1 <= 0.0498);*-;
    xl1cat_5_2 = (0.0498 < explag1 <= 0.2231);*-;
    xl1cat_5_3 = (0.2231 < explag1 <= 0.6065);*-;
    xl1cat_5_4 = (0.6065 < explag1 <= 8000);*-;
    *lagged exposure categories;
    sil2_4_1 = (0.0291 < cumexplag2 <= 1.2647);*25%-50%;
    sil2_4_2 = (1.2647 < cumexplag2 <= 6.1493);*50%-75%;
    sil2_4_3 = (6.1493 < cumexplag2 <= 235000);*75%-100%;
    * exposure window;
    exp_1_15 = cumexplag1 - cumexplag15;
    exp_2_15 = cumexplag2 - cumexplag15; *lagged additional year for exposure prediction;
    * employment, years employed/exposed (note that 'exposed' is used loosely here and considered same as 'employed' even though there are unexposed, but employed person years);
    atwork = (yrsexp>0);
    atworklag1 = (yrsexplag1>0);  
    cyxl1sq = cumyrsexplag1*cumyrsexplag1;
    cyxl1cu = cumyrsexplag1*cumyrsexplag1*cumyrsexplag1;
    cye1_3_1 = (7.3808 < cumyrsexplag1 <= 18.9727);*33.33333%-66.66667%;
    cye1_3_2 = (18.9727 < cumyrsexplag1 <= 68000);*66.66667%-100%;
%MEND;


