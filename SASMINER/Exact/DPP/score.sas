%macro score;
/* Delete Code Modifying Exported Metadata */
filename tempd "&EM_FILE_CDELTA_TRAIN";
data _null_;
if fexist('tempd') then
rc=fdelete('tempd');
run;
/*%if (%upcase("&EM_PROPERTY_METHOD") ne "NONE") and*/
/*(%upcase("&EM_PROPERTY_EXCLUDEDVARIABLE") ne "NONE")*/
/*%then %do;*/
/*%em_getname(key=OUTEST, type=DATA);*/
/*proc transpose data=&EM_USER_OUTEST*/
/*out=temp(where=(Col1 eq .));*/
/*run;*/
/*data _null_;*/
/*file tempd;*/
/*length String $200;*/
/*set temp end=eof;*/
/*if _N_=1 then put 'if upcase(NAME) in(';*/
/*string = quote(strip(upcase(_NAME_)));*/
/*put string;*/
/*if eof then do;*/
/*%if %upcase("&EM_PROPERTY_EXCLUDEDVARIABLE") eq "REJECT"*/
/*%then %do;*/
/*put ') then ROLE="REJECTED";';*/
/*%end;*/
/*%else %do;*/
/*put ') then delete;';*/
/*%end;*/
/*end;*/
/*run;*/
/*%end;*/
/*filename tempd;*/
%mend score;