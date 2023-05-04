DATA_SECTION
  init_int styr
  init_int endyr
	!! cout<<endyr<<endl;
  init_int recage
  init_int trmage
  vector age_vector(recage,trmage)
  !! for (int j=recage;j<=trmage;j++) age_vector(j) = double(j);

  init_number sigmar_in 
  number sigmarsq_in 
  !!     sigmarsq_in = sigmar_in * sigmar_in;

  init_number F_prior_mu
  init_number F_prior_sd
  number F_prior_var 
  !!     F_prior_var = F_prior_sd * F_prior_sd;
  init_number F_regularity_sd
  number F_regularity_pen
  !!     F_regularity_pen = .5/(F_regularity_sd * F_regularity_sd);

  init_int phase_M
  init_number M_prior_mu
  init_number M_prior_sd
  number M_prior_var 
  !! M_prior_var = M_prior_sd * M_prior_sd;

  init_int phase_sel
  int phase_sel2
  !!if (phase_sel>0) phase_sel2 = phase_sel+1; else phase_sel2=phase_sel;
  init_number slp_prior_mu
  init_number slp_prior_sd
  number slp_prior_var
  !! slp_prior_var = slp_prior_sd * slp_prior_sd;
  init_number a50_prior_mu
  init_number a50_prior_sd
  number a50_prior_var
  !! a50_prior_var = a50_prior_sd * a50_prior_sd;

  init_int phase_q
  init_number q_prior_mu
  init_number q_prior_sd
  number q_prior_var
  !! q_prior_var = q_prior_sd * q_prior_sd;

  init_int phase_q_rw
  init_number q_rw_sd
  number q_rw_var
  !! q_rw_var = q_rw_sd * q_rw_sd;

  init_int n_fut_yrs
  init_number fut_expl_M
  int styr_fut
  int endyr_fut
  !!styr_fut=endyr + 1;
  !!endyr_fut=styr_fut + n_fut_yrs;


  init_vector wt(recage,trmage)
  init_vector maturity(recage,trmage)

  init_vector obs_catch(styr,endyr)
  init_vector obs_catch_cv(styr,endyr)

  init_int nyrs_c_age
	number offset
	!! cout << obs_catch_cv<<endl<<nyrs_c_age<<endl;//exit(1);
  init_ivector yrs_c_age(1,nyrs_c_age)
	!! cout << yrs_c_age<<endl<<nyrs_c_age<<endl;//exit(1);
  init_matrix obs_c_age(1,nyrs_c_age,recage,trmage)
  !! if (nyrs_c_age>0) cout << "Last year of catch-age data: " <<endl<<yrs_c_age(nyrs_c_age)<<" ";
  !! if (nyrs_c_age>0) cout << obs_c_age(nyrs_c_age) <<endl;  // Print last row to screen for checking...
  // !! int xx; cout << " hit key and enter to continue" <<endl; cin>>xx;

  init_int nyrs_trend
  init_number trend_mo
  number trend_seas
  !! trend_seas = (trend_mo -1)/12;
  init_ivector yrs_trend(1,nyrs_trend)
  init_vector obs_trend(1,nyrs_trend)
  init_vector obs_trend_cv(1,nyrs_trend)

INITIALIZATION_SECTION
  F_annual .01
  catchability q_prior_mu;
  lnR 9.
  M M_prior_mu;
  sel_a50 4.
  sel_slp .4
  
PARAMETER_SECTION
  init_number lnR(1) 
  init_bounded_number sel_a50(0.,6.,phase_sel) 
  init_bounded_number  sel_slp(0.,2.,phase_sel2) 
  init_bounded_number M(0.01,1,phase_M) 
  init_bounded_vector F_annual(styr,endyr,0,1.,2);

  init_bounded_dev_vector init_dev(recage+1,trmage,-10,10,2);

  init_bounded_number catchability(0.0001,3,phase_q) 
  init_bounded_dev_vector q_rw(1,nyrs_trend,-5,5,phase_q_rw) 
  init_bounded_dev_vector rec_dev(styr,endyr,-10,10,2);
  init_bounded_dev_vector rec_dev_fut(styr_fut,endyr_fut,-10,10,5);
  number sigmarsq;
  // sdreport_number sigmar
  number sigmar

  vector sel(recage,trmage)
  vector pred_catch(styr,endyr)
  matrix pred_c_age(1,nyrs_c_age,recage,trmage)

  matrix N(styr,endyr,recage,trmage)
  matrix C(styr,endyr,recage,trmage)
  matrix F(styr,endyr,recage,trmage)
  matrix Z(styr,endyr,recage,trmage)
  matrix S(styr,endyr,recage,trmage)

  matrix N_fut(styr_fut,endyr_fut,recage,trmage)
  matrix C_fut(styr_fut,endyr_fut,recage,trmage)
  vector pred_catch_fut(styr_fut,endyr_fut)

  vector pred_trend(1,nyrs_trend)

  sdreport_vector B(styr,endyr)
  sdreport_vector Fmort(styr,endyr)
  sdreport_vector SSB(styr,endyr)
  sdreport_vector rec(styr,endyr)
  sdreport_vector B_fut(styr_fut,endyr_fut)

  vector prior(1,5)
  vector like(1,5)
  objective_function_value obj_fun

PROCEDURE_SECTION
  Compute_Mortality();
  Compute_N();
  Compute_Obj_Fun();
  if (last_phase())
	{
		do_output_vars();
    Future_projections();
	}

FUNCTION Compute_Mortality
 // Selectivity 
  sel = 1./(1. + mfexp(sel_slp*(sel_a50 - age_vector)));

 // Set mortality 
  int i;
  for (i=styr;i<=endyr;i++)
  {
     F(i) = F_annual(i) * sel ;
  }
  Z = F + M ;
  S = mfexp(-Z);
  
FUNCTION Compute_N
  // Do numbers-at-age
  //   first Year

  int i,j;
  N(styr,recage) = mfexp(lnR + rec_dev(styr));
  for ( j = recage+1 ;j <= trmage;j++)
    N(styr,j) = N(styr,j-1) * mfexp(-M+init_dev(j));

  // Plus group in first year
  N(styr,trmage) /= (1-mfexp(-M));

  //   subsequent years
  for ( i=styr+1 ; i<=endyr ; i++)
    N(i,recage) = exp(lnR + rec_dev(i));
  
  for ( i=styr+1 ; i<=endyr ; i++)
  {
    for ( j = recage+1 ;j <= trmage;j++)
    {
      N(i,j) = N(i-1,j-1) * S(i-1,j-1);
    }
    N(i,trmage) += N(i-1,trmage) * S(i-1,trmage);
  }

  // Compute Catch
  C = elem_div( elem_prod( elem_prod(N , F), 1.-S), Z);

  // Compute predicted values
  for ( i=styr ; i<=endyr ; i++)
  {
    pred_catch(i) = C(i) * wt ; 
  }
  if (nyrs_trend>0)
  {
    pred_trend.initialize();
    for ( i=1 ; i<=nyrs_trend ; i++)
      pred_trend(i) = catchability * mfexp(q_rw(i)) * elem_prod(N(yrs_trend(i)) , pow(S(yrs_trend(i)),trend_seas)) * elem_prod(wt , sel) ; // Scalar = vector * vector ;
  }

  if (nyrs_c_age>0)
    for ( i=1 ; i<=nyrs_c_age ; i++)
      pred_c_age(i) = C(yrs_c_age(i)); 

FUNCTION Compute_Obj_Fun
  // Compute objective function...
  like.initialize();
  like(1) = norm2((log(pred_catch + 0.001)-log(obs_catch +0.001)))*200; // Catch Biomass part

  if (nyrs_c_age>0)
    like(2) = 0.5*size_count(pred_c_age) * log(sum(elem_div(square(obs_c_age-pred_c_age),pred_c_age))) ;              // Catch-age likelihood 

  if (nyrs_trend>0)
    like(3) = norm2(elem_div(log(obs_trend+0.001)-log(pred_trend+0.001), obs_trend_cv));  // likelihood for the vulnerable biomass

  if (active(rec_dev))
    like(4) = norm2(rec_dev)/(2*sigmarsq_in);  
  like(4)  += norm2(init_dev);

  if (last_phase())  
  {
    // Future recruitment variability equal to historical...
    sigmarsq = sigmarsq_in ; // norm2(rec_dev)/size_count(rec_dev); 
    sigmar   = sqrt(sigmarsq); 
    like(4) += .5* norm2(rec_dev_fut)/sigmarsq ;
  }

  prior.initialize();
  prior(1)  = norm2(log(F_annual)-log(F_prior_mu)) / (2*F_prior_var);     
  if(last_phase()) 
    prior(1) += norm2(log(F)-log(mean(F))) * F_regularity_pen;     
  else
    prior(1) += .1* norm2(log(F)-log(mean(F))) * F_regularity_pen;     

  if (active(M))
    prior(2)  = square(log(M/M_prior_mu)) / (2*M_prior_var);     
  if (active(sel_slp))
  {
    prior(3)  = square(log(sel_slp/slp_prior_mu)) / (2*slp_prior_var);    
    prior(3) += square(log(sel_a50/a50_prior_mu)) / (2*a50_prior_var);     
  }
  if (active(catchability))
    prior(4)  = square(log(catchability/q_prior_mu)) / (2*q_prior_var);    
  if (active(q_rw))
    prior(5)  = norm2(first_difference(q_rw)) / (2*q_rw_var);     

  obj_fun  = sum(like);
  obj_fun += sum(prior);


FUNCTION do_output_vars
  for ( int i=styr ; i<=endyr ; i++)
	{
    rec(i)        = N(i,1) ; 
    B(i)          = N(i) * wt ; 
    SSB(i)        = N(i) * elem_prod(maturity,wt)/2 ; 
    Fmort(i)      = mean(F(i)(3,8));
	}

FUNCTION Future_projections
  int i,j;
  for ( j = recage+1 ;j <= trmage;j++)
    N_fut(styr_fut,j) = N(endyr,j-1) * mfexp(-(M+F(endyr,j-1)));

  dvar_vector F_fut = fut_expl_M * M * sel; 
  for (i=styr_fut ; i<=endyr_fut ; i++)
     N_fut(i,recage) = exp(lnR + rec_dev_fut(i));
 
  for ( i=styr_fut+1 ; i<=endyr_fut ; i++)
  {
    for ( j = recage+1 ;j <= trmage;j++)
      N_fut(i,j) = N_fut(i-1,j-1) * mfexp(-(M+F_fut(j-1))); 
 
    N_fut(i,trmage) += N_fut(i-1,trmage) * mfexp(-(M+F_fut(trmage)));
  }

  // Compute Catch
  for ( i=styr_fut ; i<=endyr_fut ; i++)
  {
    C_fut(i) = elem_prod(elem_prod(N_fut(i) , F_fut) , elem_div((1-exp(-(M+F_fut))) , (M + F_fut)));
    B_fut(i) = N_fut(i) * wt;
    pred_catch_fut(i) = C_fut(i) * wt;
  }


FUNCTION write_R
  REPORT(like);
  REPORT(prior);
  REPORT(B);
  REPORT(F_annual);
  REPORT(sel);
  REPORT(obs_c_age);
  REPORT(pred_c_age);
  REPORT(N);
  REPORT(N_fut);
  REPORT(yrs_trend);
  REPORT(obs_trend);
  REPORT(pred_trend);
  dvar_vector Q = catchability * mfexp(q_rw);
  REPORT(Q);
  REPORT(B_fut)       ;
  REPORT(SSB)       ;
  REPORT(Fmort)       ;
  REPORT(rec)       ;
  REPORT(pred_catch_fut);

REPORT_SECTION
  save_gradients(gradients);
  if (last_phase())
    write_R();
  cout <<" Just finished phase " <<current_phase()<<" +++++++++++++++++++++++++++++++++++"<<endl;
  report<<"Stock summary   "<<endl;
  report<<endl<<"Selectivity "<<endl<< age_vector<<endl<<sel<<endl<<endl;
  report<<"Year Biomass F Catch/Biomass SpBiom Recruit"<<endl;
  for (int i=styr;i<=endyr;i++) 
    report << i                 <<" "
           <<B(i)               <<" "
           <<F_annual(i)        <<" "
           <<pred_catch(i)/B(i) <<" "
           <<N(i)*elem_prod(maturity,wt)<<" "
           <<N(i,recage)        <<" "
           <<endl;
  report<<"Numbers_at_age"<<endl;
  report<<"Year "<<age_vector<<endl; for (int i=styr;i<=endyr;i++) report << i <<" " <<N(i) <<endl;
  report<<"Catch_at_age"<<endl;
  report<<"Year "<<age_vector<<endl; for (int i=styr;i<=endyr;i++) report << i <<" " <<C(i) <<endl;
  report<<"F_at_age"<<endl;
  report<<"Year "<<age_vector<<endl; for (int i=styr;i<=endyr;i++) report << i <<" " <<F(i) <<endl;

  report<<endl<<"Likelihoods"<<endl<< 
          " Catch catch_age survey recruitment" <<endl<<
          like<<endl;
  report<<endl<<"Priors"<<endl<< 
          " F M Sel q q_trend"<< endl<<
            prior             << endl;
  report<<endl<<"Trend/Index data"<<endl;
  if (nyrs_trend>0)
  {
    report<<endl<< "Fits to data "<<endl;
    report<<"Year Obs Predicted catchability"<<endl;
    for (int i=1;i<=nyrs_trend;i++) 
      report << yrs_trend(i)<<" "<<obs_trend(i)<<" "
             <<pred_trend(i)<<" "
             <<catchability * mfexp(q_rw(i)) <<" "<<
             endl;
  }
  report<<endl<<"Catch Biomass "<<endl<< "Obs "<<obs_catch<<endl<<"Pred "<<pred_catch<<endl;
  if (nyrs_c_age>0)
  {
    report<<endl<<"Catch-at-age"<<endl<<"Observed "<<endl<< "Year "<<age_vector<<endl;
    for (int i=1;i<=nyrs_c_age;i++) report<<yrs_c_age(i)<<"  "<<obs_c_age(i)<<endl;
    report<<endl<<"Catch-at-age"<<endl<<"Predicted "<<endl<< "Year "<<age_vector<<endl;
    for (int i=1;i<=nyrs_c_age;i++) report<<yrs_c_age(i)<<"  "<<pred_c_age(i)<<endl;
  }
  
  report<<"Future "<<styr_fut<<" - "<<endyr_fut<<endl;
  report<<"Year Biomass Catch SpBiom Recruit"<<endl;
  for (int i=styr_fut;i<=endyr_fut;i++) 
    report << i                     <<" "
           <<B_fut(i)               <<" "
           <<pred_catch_fut(i)      <<" "
           <<N_fut(i)*elem_prod(maturity,wt)<<" "
           <<endl;
  report<<"Future_N_at_age"<<endl;
  report<<"Year "<<age_vector<<endl; for (int i=styr_fut;i<=endyr_fut;i++) report << i <<" " <<N_fut(i) <<endl;
  report<<"Future_C_at_age"<<endl;
  report<<"Year "<<age_vector<<endl; for (int i=styr_fut;i<=endyr_fut;i++) report << i <<" " <<C_fut(i) <<endl;

FINAL_SECTION
  REPORT(Fmort.sd);
  REPORT(SSB.sd);
  REPORT(rec.sd);

GLOBALS_SECTION
  #include <admodel.h>
  #include <time.h>
  // Define objects for report file, echoinput, etc.
  /**
  \def report(object)
  Prints name and value of \a object on ADMB report %ofstream file.
  */
  #undef REPORT
  #define REPORT(object) Rreport <<  #object "\n" << setw(8)  << setprecision(4) << setfixed() << object << endl;
  ofstream Rreport("R_report.rep");
