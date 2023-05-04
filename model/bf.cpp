#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
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
#ifdef DEBUG
  #include <chrono>
#endif
#include <admodel.h>
#ifdef USE_ADMB_CONTRIBS
#include <contrib.h>

#endif
  extern "C"  {
    void ad_boundf(int i);
  }
#include <bf.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  adstring tmpstring;
  tmpstring=adprogram_name + adstring(".dat");
  if (argc > 1)
  {
    int on=0;
    if ( (on=option_match(argc,argv,"-ind"))>-1)
    {
      if (on>argc-2 || argv[on+1][0] == '-')
      {
        cerr << "Invalid input data command line option"
                " -- ignored" << endl;
      }
      else
      {
        tmpstring = adstring(argv[on+1]);
      }
    }
  }
  global_datafile = new cifstream(tmpstring);
  if (!global_datafile)
  {
    cerr << "Error: Unable to allocate global_datafile in model_data constructor.";
    ad_exit(1);
  }
  if (!(*global_datafile))
  {
    delete global_datafile;
    global_datafile=NULL;
  }
  styr.allocate("styr");
  endyr.allocate("endyr");
 cout<<endyr<<endl;
  recage.allocate("recage");
  trmage.allocate("trmage");
  age_vector.allocate(recage,trmage);
 for (int j=recage;j<=trmage;j++) age_vector(j) = double(j);
  sigmar_in.allocate("sigmar_in");
     sigmarsq_in = sigmar_in * sigmar_in;
  F_prior_mu.allocate("F_prior_mu");
  F_prior_sd.allocate("F_prior_sd");
     F_prior_var = F_prior_sd * F_prior_sd;
  F_regularity_sd.allocate("F_regularity_sd");
     F_regularity_pen = .5/(F_regularity_sd * F_regularity_sd);
  phase_M.allocate("phase_M");
  M_prior_mu.allocate("M_prior_mu");
  M_prior_sd.allocate("M_prior_sd");
 M_prior_var = M_prior_sd * M_prior_sd;
  phase_sel.allocate("phase_sel");
if (phase_sel>0) phase_sel2 = phase_sel+1; else phase_sel2=phase_sel;
  slp_prior_mu.allocate("slp_prior_mu");
  slp_prior_sd.allocate("slp_prior_sd");
 slp_prior_var = slp_prior_sd * slp_prior_sd;
  a50_prior_mu.allocate("a50_prior_mu");
  a50_prior_sd.allocate("a50_prior_sd");
 a50_prior_var = a50_prior_sd * a50_prior_sd;
  phase_q.allocate("phase_q");
  q_prior_mu.allocate("q_prior_mu");
  q_prior_sd.allocate("q_prior_sd");
 q_prior_var = q_prior_sd * q_prior_sd;
  phase_q_rw.allocate("phase_q_rw");
  q_rw_sd.allocate("q_rw_sd");
 q_rw_var = q_rw_sd * q_rw_sd;
  n_fut_yrs.allocate("n_fut_yrs");
  fut_expl_M.allocate("fut_expl_M");
styr_fut=endyr + 1;
endyr_fut=styr_fut + n_fut_yrs;
  wt.allocate(recage,trmage,"wt");
  maturity.allocate(recage,trmage,"maturity");
  obs_catch.allocate(styr,endyr,"obs_catch");
  obs_catch_cv.allocate(styr,endyr,"obs_catch_cv");
  nyrs_c_age.allocate("nyrs_c_age");
 cout << obs_catch_cv<<endl<<nyrs_c_age<<endl;//exit(1);
  yrs_c_age.allocate(1,nyrs_c_age,"yrs_c_age");
 cout << yrs_c_age<<endl<<nyrs_c_age<<endl;//exit(1);
  obs_c_age.allocate(1,nyrs_c_age,recage,trmage,"obs_c_age");
 if (nyrs_c_age>0) cout << "Last year of catch-age data: " <<endl<<yrs_c_age(nyrs_c_age)<<" ";
 if (nyrs_c_age>0) cout << obs_c_age(nyrs_c_age) <<endl;  // Print last row to screen for checking...
  nyrs_trend.allocate("nyrs_trend");
  trend_mo.allocate("trend_mo");
 trend_seas = (trend_mo -1)/12;
  yrs_trend.allocate(1,nyrs_trend,"yrs_trend");
  obs_trend.allocate(1,nyrs_trend,"obs_trend");
  obs_trend_cv.allocate(1,nyrs_trend,"obs_trend_cv");
}

void model_parameters::initializationfunction(void)
{
  F_annual.set_initial_value(.01);
  catchability.set_initial_value(q_prior_mu);
  lnR.set_initial_value(9.);
  M.set_initial_value(M_prior_mu);
  sel_a50.set_initial_value(4.);
  sel_slp.set_initial_value(.4);
  if (global_datafile)
  {
    delete global_datafile;
    global_datafile = NULL;
  }
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  lnR.allocate(1,"lnR");
  sel_a50.allocate(0.,6.,phase_sel,"sel_a50");
  sel_slp.allocate(0.,2.,phase_sel2,"sel_slp");
  M.allocate(0.01,1,phase_M,"M");
  F_annual.allocate(styr,endyr,0,1.,2,"F_annual");
  init_dev.allocate(recage+1,trmage,-10,10,2,"init_dev");
  catchability.allocate(0.0001,3,phase_q,"catchability");
  q_rw.allocate(1,nyrs_trend,-5,5,phase_q_rw,"q_rw");
  rec_dev.allocate(styr,endyr,-10,10,2,"rec_dev");
  rec_dev_fut.allocate(styr_fut,endyr_fut,-10,10,5,"rec_dev_fut");
  sigmarsq.allocate("sigmarsq");
  #ifndef NO_AD_INITIALIZE
  sigmarsq.initialize();
  #endif
  sigmar.allocate("sigmar");
  #ifndef NO_AD_INITIALIZE
  sigmar.initialize();
  #endif
  sel.allocate(recage,trmage,"sel");
  #ifndef NO_AD_INITIALIZE
    sel.initialize();
  #endif
  pred_catch.allocate(styr,endyr,"pred_catch");
  #ifndef NO_AD_INITIALIZE
    pred_catch.initialize();
  #endif
  pred_c_age.allocate(1,nyrs_c_age,recage,trmage,"pred_c_age");
  #ifndef NO_AD_INITIALIZE
    pred_c_age.initialize();
  #endif
  N.allocate(styr,endyr,recage,trmage,"N");
  #ifndef NO_AD_INITIALIZE
    N.initialize();
  #endif
  C.allocate(styr,endyr,recage,trmage,"C");
  #ifndef NO_AD_INITIALIZE
    C.initialize();
  #endif
  F.allocate(styr,endyr,recage,trmage,"F");
  #ifndef NO_AD_INITIALIZE
    F.initialize();
  #endif
  Z.allocate(styr,endyr,recage,trmage,"Z");
  #ifndef NO_AD_INITIALIZE
    Z.initialize();
  #endif
  S.allocate(styr,endyr,recage,trmage,"S");
  #ifndef NO_AD_INITIALIZE
    S.initialize();
  #endif
  N_fut.allocate(styr_fut,endyr_fut,recage,trmage,"N_fut");
  #ifndef NO_AD_INITIALIZE
    N_fut.initialize();
  #endif
  C_fut.allocate(styr_fut,endyr_fut,recage,trmage,"C_fut");
  #ifndef NO_AD_INITIALIZE
    C_fut.initialize();
  #endif
  pred_catch_fut.allocate(styr_fut,endyr_fut,"pred_catch_fut");
  #ifndef NO_AD_INITIALIZE
    pred_catch_fut.initialize();
  #endif
  pred_trend.allocate(1,nyrs_trend,"pred_trend");
  #ifndef NO_AD_INITIALIZE
    pred_trend.initialize();
  #endif
  B.allocate(styr,endyr,"B");
  Fmort.allocate(styr,endyr,"Fmort");
  SSB.allocate(styr,endyr,"SSB");
  rec.allocate(styr,endyr,"rec");
  B_fut.allocate(styr_fut,endyr_fut,"B_fut");
  prior.allocate(1,5,"prior");
  #ifndef NO_AD_INITIALIZE
    prior.initialize();
  #endif
  like.allocate(1,5,"like");
  #ifndef NO_AD_INITIALIZE
    like.initialize();
  #endif
  obj_fun.allocate("obj_fun");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
}

void model_parameters::userfunction(void)
{
  obj_fun =0.0;
  Compute_Mortality();
  Compute_N();
  Compute_Obj_Fun();
  if (last_phase())
	{
		do_output_vars();
    Future_projections();
	}
#ifdef DEBUG
  std::cout << "DEBUG: Total gradient stack used is " << gradient_structure::get()->GRAD_STACK1->total() << " out of " << gradient_structure::get_GRADSTACK_BUFFER_SIZE() << std::endl;;
  std::cout << "DEBUG: Total dvariable address used is " << gradient_structure::get()->GRAD_LIST->total_addresses() << " out of " << gradient_structure::get_MAX_DLINKS() << std::endl;;
  std::cout << "DEBUG: Total dvariable address used is " << gradient_structure::get()->ARR_LIST1->get_max_last_offset() << " out of " << gradient_structure::get_ARRAY_MEMBLOCK_SIZE() << std::endl;;
#endif
}

void model_parameters::Compute_Mortality(void)
{
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
  
}

void model_parameters::Compute_N(void)
{
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
}

void model_parameters::Compute_Obj_Fun(void)
{
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
}

void model_parameters::do_output_vars(void)
{
  for ( int i=styr ; i<=endyr ; i++)
	{
    rec(i)        = N(i,1) ; 
    B(i)          = N(i) * wt ; 
    SSB(i)        = N(i) * elem_prod(maturity,wt)/2 ; 
    Fmort(i)      = mean(F(i)(3,8));
	}
}

void model_parameters::Future_projections(void)
{
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
}

void model_parameters::write_R(void)
{
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
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
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
}

void model_parameters::final_calcs()
{
  REPORT(Fmort.sd);
  REPORT(SSB.sd);
  REPORT(rec.sd);
}

void model_parameters::preliminary_calculations(void){
#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::set_runtime(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
    gradient_structure::set_NO_DERIVATIVES();
#ifdef DEBUG
  #ifndef __SUNPRO_C
std::feclearexcept(FE_ALL_EXCEPT);
  #endif
  auto start = std::chrono::high_resolution_clock::now();
#endif
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint = defaults::iprint;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
#ifdef DEBUG
  std::cout << endl << argv[0] << " elapsed time is " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count() << " microseconds." << endl;
  #ifndef __SUNPRO_C
bool failedtest = false;
if (std::fetestexcept(FE_DIVBYZERO))
  { cerr << "Error: Detected division by zero." << endl; failedtest = true; }
if (std::fetestexcept(FE_INVALID))
  { cerr << "Error: Detected invalid argument." << endl; failedtest = true; }
if (std::fetestexcept(FE_OVERFLOW))
  { cerr << "Error: Detected overflow." << endl; failedtest = true; }
if (std::fetestexcept(FE_UNDERFLOW))
  { cerr << "Error: Detected underflow." << endl; }
if (failedtest) { std::abort(); } 
  #endif
#endif
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
