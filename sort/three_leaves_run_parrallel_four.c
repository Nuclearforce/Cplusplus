#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
#include <TH1.h>
#include <TH2.h>
#include <TROOT.h>
#include <iostream>
#include <vector>
#include <TRandom3.h>
#include <cstdio>
#include "TSystem.h"

//start 'global'
TRandom3 *r = new TRandom3();
TFile f_hist("../histogram/14C/04_25_2018_82.root", "recreate");
//TFile f_hist("../histogram/13C_uncalibrated_before_experiment.root", "recreate");

TChain DataChain("data");
TTreeReader myReader(&DataChain);
TTreeReaderArray<unsigned int> TRenergy(myReader, "energy");
TTreeReaderArray<unsigned int> TRtime(myReader, "time");
TTreeReaderArray<unsigned int> TRmodchan(myReader, "modchan");

TChain Background("data");
TTreeReader myReader_background(&Background);
TTreeReaderArray<unsigned int> background_energy(myReader_background, "energy");
TTreeReaderArray<unsigned int> background_time(myReader_background, "time");
TTreeReaderArray<unsigned int> background_modchan(myReader_background, "modchan");

//extra
Bool_t for_testing = false;

//gate
Bool_t gate_on_excitation_spectra;
Bool_t particle_gamma;
Bool_t energy_sharing_ring;
Bool_t energy_sharing_sector;
Bool_t ring_sector_time_gate;
Double_t ring_time;
Double_t sector_time;
Double_t gamma_time;
Bool_t bad_gamma;

//coincidence gate
/*Double_t lower_coinc = 912;
 Double_t upper_coinc = 922;*/
Int_t print_clover_count = 0;
Int_t print_oslo_count = 0;
Int_t print_army_count = 0;

//lower and upper limit for excitation spectrum gate
Double_t lower_lim = 40000;
Double_t upper_lim = 46000;

Int_t ring_counter=0;
Int_t sector_counter=0;

//histograms
TH1F* hPixieADC[4][16];
//TH1F* hPixieADC_eff[4][16];
//TH1F* time_background_eff[4][16];
TH1F* hPixieADC_ab[4][16];
TH1F* hPixieADC_subtracted[4][16];
TH1F* hPixieADC_pg[4][16];
TH1F* backgroundADC[4][16];
TH1F* time_background_not_dopler[4][16];
TH1F* time_background_not_dopler_pg[4][16];
TH1F* time_Ge;
TH1F* time_Cl;
TH1F* time_La;
TH1F* time_La_army;
TH1F* pg_time;
TH1F* ring_sector_time;
TH1F* time_background[4][16];
TH1F* time_background_pg[4][16];
TH1F* time_Si;
TH1F* Clover1;
TH1F* Clover3;
TH1F* Clover1ab;
TH1F* Clover3ab;
TH1F* rings_hit;
TH1F* all_detectors_complete;
TH1F* all_detectors_ab;
TH1F* all_Ge_ab;
TH1F* all_LaBr_ab;
TH1F* all_Ge;
TH1F* all_LaBr;
TH1F* one_ring_singles;
TH1F* one_ring_timegated;
TH1F* one_ring_timegated_pggated;
TH2F *ring_en_sector_en;
TH2F *ring_en_gamma_en;
TH2F *ring_en_gamma_en_background;
TH1I* rings_fired;
TH1I* sectors_fired;
TH1F* ring_middel;
TH1F* ring_up;
TH1F* ring_down;

//calibration
Double_t calib_a[4][16];
Double_t calib_b[4][16];
Double_t calib_c[4][16];
Double_t calib_d[4][16];
Double_t calib_e[4][16];

//Doppler correction
Double_t beta = 0.0;
Double_t theta[4][16];
Double_t phi[4][16];
Int_t ring_mod;
Int_t ring_chn;
Int_t sector_mod;
Int_t sector_chn;

//beta calculation
Double_t beta_ring_en=0.0;
Double_t ring_en_for_hist=0.0;
Double_t sector_en_for_hist=0.0;

//energy_sharing
Double_t energy_sharing_ring_en=0;
Double_t energy_sharing_ring_chn_use=0;
Double_t energy_sharing_ring_mod_use=0;
Double_t energy_sharing_ring_chn_ignore=0;
Double_t energy_sharing_ring_mod_ignore=0;
Double_t energy_sharing_sector_en=0;
Double_t energy_sharing_sector_chn_use=0;
Double_t energy_sharing_sector_chn_ignore=0;
Double_t energy_share_use_time_ring=99999;
Double_t energy_share_use_time_sector=99999;
Bool_t bad_particle;

//Add-back
Double_t ab_use_energy_cl1;
Double_t ab_use_crystal_cl1;
Double_t ab_ignore_crystal_cl1;
Bool_t ab_event_cl1=false;
Double_t ab_use_energy_cl3;
Double_t ab_use_crystal_cl3;
Double_t ab_ignore_crystal_cl3;
Bool_t ab_event_cl3=false;
Double_t test_counter=0;
Double_t add_back_use_time1=99999;
Double_t add_back_use_time3=99999;

//coincidence
Int_t coincidence_counter=0;
Double_t gamma_2_time = 0;

//function declaration
void init_histos();
void read_calibrate_file();
void make_chain();
void first_loop();
void second_loop();
Double_t do_calibrate(Double_t raw_energy, Int_t module, Int_t channel);
Double_t do_doppler(Double_t en,Int_t gamma_mod, Int_t gamma_chn);
Double_t do_beta(Double_t ring_en);

void three_leaves_run_parrallel_four()
{
    
    make_chain();//add trees to chain
    init_histos(); //initialise histograms
    read_calibrate_file(); //read calibration file
    
    // start timer for progress bar
    std::clock_t start;
    Double_t duration;
    Long64_t timer = 0;
    Long64_t max_event = DataChain.GetEntries();
    
    Bool_t stop_now = false;
    Long64_t max_event_testing = 1e7;
    start = std::clock();
    
    stop_now = false;
    
    ////////////////////// Loop over all entries of the TChain.////////////////////////////
    while (myReader.Next()) {
        first_loop();
        //cout<<"one"<<endl;
        //gate_on_excitation_spectra=true; //set to true if no excitation gate to be used
        //bad_particle=false;//set to false if no Si data
        //bad_gamma=false;//set to false if no Si data
        
        /* if(bad_gamma==true||gate_on_excitation_spectra==false){ ///////////////////used for testing
         //skip event
         }*/
        
        if(bad_particle==true||bad_gamma==true||gate_on_excitation_spectra==false){
            //skip event
        }
        
        else if (gate_on_excitation_spectra==true){
            //cout<<"good event" <<endl;
            //cout<<"spam"<<endl;
            second_loop();
            //cout<<"two"<<endl;
        }
        
        
        if (for_testing==true){
            if(timer%5000000==0){
                cout << "Event #: " << timer  << " Progress: " << 100.0*timer/max_event_testing << " [%] \r";
                cout.flush();
            }
            timer++;
            if(timer==max_event_testing){
                stop_now=true;
            }
        }
        else{
            if(timer%50000000==0){
                cout << "Event #: " << timer  << " Progress: " << 100.0*timer/max_event << " [%] \r";
                cout.flush();
            }
            timer++;
        }
        if (stop_now==true){
            break;
        }
    }
    
    ////////////////write and close histogram file////////////////////
    Int_t j=0;
    Int_t i=0;
    Double_t Labr_subtract=13.0;
    Double_t Ge_subtract=12.0;
    Double_t background_size=800.0;
    for(j=0;j<4;j++){
        for(i=0;i<16;i++){
            if (j==2&&(i==0||i==1||i==2||i==3||i==4||i==5)){
                hPixieADC_pg[j][i]->Add(time_background[j][i],-(Labr_subtract/background_size));
                hPixieADC_ab[j][i]->Add(time_background_not_dopler[j][i],-(Labr_subtract/background_size));
                hPixieADC_pg[j][i]->Add(time_background_pg[j][i],-(Labr_subtract/background_size));
                hPixieADC_ab[j][i]->Add(time_background_not_dopler_pg[j][i],-(Labr_subtract/background_size));
                //hPixieADC_eff[j][i]->Add(time_background_eff[j][i],-(16.0/600.0));
            }
            else{
                hPixieADC_pg[j][i]->Add(time_background[j][i],-(Ge_subtract/background_size));
                hPixieADC_ab[j][i]->Add(time_background_not_dopler[j][i],-(Ge_subtract/background_size));
                hPixieADC_pg[j][i]->Add(time_background_pg[j][i],-(Ge_subtract/background_size));
                hPixieADC_ab[j][i]->Add(time_background_not_dopler_pg[j][i],-(Ge_subtract/background_size));
                //hPixieADC_eff[j][i]->Add(time_background_eff[j][i],-(16.0/600.0));
            }
        }
    }
    Clover1->Add(hPixieADC_pg[1][0],1);
    Clover1->Add(hPixieADC_pg[1][1],1);
    Clover1->Add(hPixieADC_pg[1][2],1);
    Clover1->Add(hPixieADC_pg[1][3],1);
    
    Clover3->Add(hPixieADC_pg[3][0],1);
    Clover3->Add(hPixieADC_pg[3][1],1);
    Clover3->Add(hPixieADC_pg[3][2],1);
    Clover3->Add(hPixieADC_pg[3][3],1);
    
    
    Clover1ab->Add(hPixieADC_ab[1][0],1);
    Clover1ab->Add(hPixieADC_ab[1][1],1);
    Clover1ab->Add(hPixieADC_ab[1][2],1);
    Clover1ab->Add(hPixieADC_ab[1][3],1);
    
    Clover3ab->Add(hPixieADC_ab[3][0],1);
    Clover3ab->Add(hPixieADC_ab[3][1],1);
    Clover3ab->Add(hPixieADC_ab[3][2],1);
    Clover3ab->Add(hPixieADC_ab[3][3],1);
    
    all_detectors_complete->Add(Clover1,1.0);
    all_detectors_complete->Add(Clover3,1.0);
    all_detectors_complete->Add(hPixieADC_pg[2][0],1.0);
    all_detectors_complete->Add(hPixieADC_pg[2][1],1.0);
    all_detectors_complete->Add(hPixieADC_pg[2][2],1.0);
    all_detectors_complete->Add(hPixieADC_pg[2][3],1.0);
    all_detectors_complete->Add(hPixieADC_pg[3][5],1.0);
    
    all_detectors_ab->Add(Clover1ab,1.0);
    all_detectors_ab->Add(Clover3ab,1.0);
    all_detectors_ab->Add(hPixieADC_ab[2][0],1.0);
    all_detectors_ab->Add(hPixieADC_ab[2][1],1.0);
    all_detectors_ab->Add(hPixieADC_ab[2][2],1.0);
    all_detectors_ab->Add(hPixieADC_ab[2][3],1.0);
    all_detectors_ab->Add(hPixieADC_ab[3][5],1.0);
    
    all_Ge_ab->Add(hPixieADC_ab[3][5],1.0);
    all_Ge_ab->Add(Clover1ab,1.0);
    all_Ge_ab->Add(Clover3ab,1.0);
    
    all_Ge->Add(hPixieADC_pg[3][5],1.0);
    all_Ge->Add(Clover1,1.0);
    all_Ge->Add(Clover3,1.0);
    
    all_LaBr_ab->Add(hPixieADC_ab[2][0],1.0);
    all_LaBr_ab->Add(hPixieADC_ab[2][1],1.0);
    all_LaBr_ab->Add(hPixieADC_ab[2][2],1.0);
    all_LaBr_ab->Add(hPixieADC_ab[2][3],1.0);
    
    all_LaBr->Add(hPixieADC_pg[2][0],1.0);
    all_LaBr->Add(hPixieADC_pg[2][1],1.0);
    all_LaBr->Add(hPixieADC_pg[2][2],1.0);
    all_LaBr->Add(hPixieADC_pg[2][3],1.0);
    
    ring_en_gamma_en->Add(ring_en_gamma_en_background,-(Labr_subtract/background_size));
    
    //cout<<gamma_amount_of_detectors_count<<endl;
    f_hist.Write();
    duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
    cout << "Timer: " << duration << " [s] to analyze " << max_event << " events." << endl;
    f_hist.Close();
    
}

////////////////////////first loop, checks for add-back and gates//////////////////
void first_loop(){
    
    Int_t entry_id;
    Int_t entry_size;
    ring_counter=0;
    sector_counter=0;
    Int_t gamma_counter=0;
    
    Double_t time;
    Double_t en;
    Double_t calib_en;
    Int_t LaBr_counter=0;
    Int_t LaBr_gamma_chn=9999;
    
    //gates
    gate_on_excitation_spectra = false;
    Double_t gamma_amount_of_detectors_check =0;
    ring_time=999999;
    sector_time=99999;
    gamma_time=99999;
    bad_gamma=false;
    Double_t temp_labr_en = 0;
    
    //energy_sharing
    Double_t old_sharing_ring_en=0;
    Double_t old_sharing_ring_chn=0;
    Double_t old_sharing_ring_mod=0;
    Double_t old_sharing_sector_en=0;
    Double_t old_sharing_sector_chn=0;
    bad_particle=false;
    energy_sharing_ring=false;
    energy_sharing_sector=false;
    ring_sector_time_gate=false;
    
    //Add-back for clover 1
    Double_t ab_oldtime_cl1=0;
    Double_t ab_olden_cl1=0;
    Double_t ab_chn_cl1=0;
    Double_t count_cl1=0;
    Double_t ab_oldlocation_cl1=0;
    ab_event_cl1=false;
    
    //Add-back for clover 3
    Double_t ab_oldtime_cl3=0;
    Double_t ab_olden_cl3=0;
    Double_t ab_chn_cl3=0;
    Double_t count_cl3=0;
    Double_t ab_oldlocation_cl3=0;
    ab_event_cl3=false;
    
    for (entry_id = 0, entry_size = TRenergy.GetSize(); entry_id < entry_size; ++entry_id) {
        
        Int_t modchan=TRmodchan[entry_id];
        Int_t mod=modchan/16;
        Int_t chn=modchan%16;
        time=TRtime[entry_id];
        en=TRenergy[entry_id];
        
        
        
        if(mod==2&&(chn==0 || chn==1 || chn==2 || chn==3 || chn==4 || chn==5)){
            //for the LaBr we want to doppler shift first and then calibrate. For the Ge we want to calibrate, Add-back then doppler shift.
            calib_en=en;
        }
        else{
            calib_en=do_calibrate(en,mod,chn);
        }
        
        ///////////////energy sharing////////////////////////
        
        if( (mod==1 && (chn==8||chn==9||chn==10||chn==11||chn==12||chn==13||chn==14||chn==15)) || (mod==2 && (chn==8||chn==9||chn==10||chn==11||chn==12||chn==13||chn==14||chn==15)) || (mod==3 && (chn==7||chn==8||chn==9||chn==10||chn==11||chn==12||chn==13||chn==14||chn==15)) ) {
            
            
            if (ring_counter==0){
                beta_ring_en=calib_en;
                ring_en_for_hist=calib_en;
                ring_time=time;
                ring_mod = mod;
                ring_chn = chn;
                old_sharing_ring_en=calib_en;
                old_sharing_ring_chn=chn;
                old_sharing_ring_mod=mod;
                energy_share_use_time_ring=time;
            }
            else{
                bad_particle=true;
            }
            if(ring_counter==1){
                //cout <<old_sharing_ring_chn-chn<<" " << old_sharing_ring_chn<< " " << chn << endl;
                if(abs(old_sharing_ring_chn-chn)==1){
                    energy_sharing_ring_en=old_sharing_sector_en+en;
                    energy_sharing_ring_chn_use=chn;
                    energy_sharing_ring_chn_ignore=old_sharing_ring_chn;
                    energy_sharing_ring_mod_use=mod;
                    energy_sharing_ring_mod_ignore=old_sharing_ring_mod;
                    energy_sharing_ring=true;
                    ring_mod = mod;
                    ring_chn = chn;
                    beta_ring_en=energy_sharing_ring_en;
                    ring_en_for_hist=energy_sharing_ring_en;
                    if(time<energy_share_use_time_ring){
                        energy_share_use_time_ring=time;
                    }
                    bad_particle=false;
                    //cout <<old_sharing_ring_chn-chn<<" " << old_sharing_ring_chn<< " " << chn << endl;
                }
                else{
                    bad_particle=true;
                    //cout <<"bad_particle"<<endl;
                }
                if(energy_sharing_ring_en<lower_lim||energy_sharing_ring_en>upper_lim){
                    bad_particle=true;
                }
            }
            
            if(ring_counter>1||en<1){
                bad_particle=true;
                //cout <<"bad_particle"<<endl;
            }
            
            
            if((energy_sharing_ring_en<lower_lim||energy_sharing_ring_en>upper_lim)&&(energy_sharing_ring==true)){
                bad_particle=true;
            }
            
            
            ring_counter++;
            
        }
        //////////excitation spectra gate (particle spectra gates)//////////////////
        if(ring_en_for_hist>lower_lim && ring_en_for_hist<upper_lim){
            gate_on_excitation_spectra=true;
        }
        if(mod==0){
            if (sector_counter==0){
                sector_en_for_hist=calib_en;
                sector_time=time;
                old_sharing_sector_en=calib_en;
                old_sharing_sector_chn=chn;
                sector_mod = mod;
                sector_chn = chn;
                energy_share_use_time_sector=time;
            }
            else{
                bad_particle=true;
            }
            if (sector_counter==1){
                if(abs(old_sharing_sector_chn-chn)==1){
                    energy_sharing_sector_en=old_sharing_sector_en+en;
                    energy_sharing_sector_chn_use=chn;
                    energy_sharing_sector_chn_ignore=old_sharing_sector_chn;
                    energy_sharing_sector=true;
                    sector_mod = mod;
                    sector_chn = chn;
                    sector_en_for_hist=energy_sharing_sector_en;
                    if(time<energy_share_use_time_sector){
                        energy_share_use_time_sector=time;
                    }
                }
                else{
                    bad_particle=true;
                }
                if(energy_sharing_sector_en<lower_lim||energy_sharing_sector_en>upper_lim){
                    bad_particle=true;
                }
            }
            
            if(sector_counter>1 || en < 1){
                bad_particle=true;
                
            }
            if( (calib_en>lower_lim && calib_en<upper_lim)){
                gate_on_excitation_spectra=true;
                
            }
            if((energy_sharing_sector_en<lower_lim||energy_sharing_sector_en>upper_lim)&&(energy_sharing_sector==true)){
                bad_particle=true;
            }
            sector_counter++;
            
        }
        //////////excitation spectra gate (particle spectra gates)//////////////////
        if(sector_en_for_hist<lower_lim || sector_en_for_hist>upper_lim){
            gate_on_excitation_spectra=true;
        }
        if(ring_en_for_hist<lower_lim || ring_en_for_hist>upper_lim){  //this is technically redundant but I feel better having two if statements doing the same thing.
            bad_particle=true;
        }
        if(sector_en_for_hist<lower_lim || sector_en_for_hist>upper_lim){ //this is technically redundant but I feel better having two if statements doing the same thing.
            bad_particle=true;
        }
        //////////////////////add-back////////////////////
        if (mod==1 && (chn==0 || chn==1 || chn==2 || chn==3) && calib_en < 10000 && calib_en>20){
            if(count_cl1==0){
                add_back_use_time1=time;
            }
            if ((time-ab_oldtime_cl1)<101.0 && count_cl1==1){
                ab_event_cl1=true;
                ab_use_energy_cl1=ab_olden_cl1+calib_en;
                ab_use_crystal_cl1=ab_oldlocation_cl1;
                ab_ignore_crystal_cl1=chn;
                if(time<add_back_use_time1){
                    add_back_use_time1=time;
                }
                
                
            }
            if(count_cl1>2||ab_olden_cl1+calib_en>10000){
                ab_event_cl1=false;
            }
            ab_oldtime_cl1=time;
            ab_olden_cl1=calib_en;
            ab_oldlocation_cl1=chn;
            count_cl1++;
        }
        
        if (mod==3 && (chn==0 || chn==1 || chn==2 || chn==3) && calib_en < 10000 && calib_en>20){
            if(count_cl3==0){
                add_back_use_time3=time;
            }
            if ((time-ab_oldtime_cl3)<101.0 && count_cl3==1){
                ab_event_cl3=true;
                ab_use_energy_cl3=ab_olden_cl3+calib_en;
                ab_use_crystal_cl3=ab_oldlocation_cl3;
                ab_ignore_crystal_cl3=chn;
                if(time<add_back_use_time3){
                    add_back_use_time3=time;
                }
            }
            if(count_cl3>2||ab_olden_cl3+calib_en>10000){
                ab_event_cl3=false;
            }
            ab_oldtime_cl3=time;
            ab_olden_cl3=calib_en;
            ab_oldlocation_cl3=chn;
            count_cl3++;
        }
        
        if( ((mod==1 && (chn==0 || chn==1 || chn==2 || chn==3)) || (mod==3 && (chn==0 || chn==1 || chn==2 || chn==3||chn==5)) || (mod==2 && (chn==0 || chn==1 || chn==2 || chn==3 || chn==4 || chn==5))) && calib_en < 10000 && calib_en>20){
            if (gamma_counter==0){
                gamma_time=time;
            }
            if (gamma_counter==1){
                gamma_2_time=time;
            }
            
            if(ab_event_cl3==true||ab_event_cl1==true){
                bad_gamma=false;
            }
            
            if((abs(gamma_amount_of_detectors_check-mod)>0 && gamma_counter>0)||gamma_counter>2){
                bad_gamma=true;
            }
            
            if(gamma_counter>1 && ab_event_cl3==false && ab_event_cl1==false){
                
                bad_gamma=true;
            }
            
            if(mod==2 && (chn==0 || chn==1 || chn==2 || chn==3)){
                LaBr_gamma_chn=chn;
                if ((LaBr_counter>0 && abs(LaBr_gamma_chn-chn)>0)||LaBr_counter>0){
                    bad_gamma=true;
                }
                LaBr_counter++;
            }
            
            
            
            gamma_counter++;
            gamma_amount_of_detectors_check=mod;
        }
        
        if( ((mod==1 && (chn==0 || chn==1 || chn==2 || chn==3)) || (mod==3 && (chn==0 || chn==1 || chn==2 || chn==3||chn==5)) || (mod==2 && (chn==0 || chn==1 || chn==2 || chn==3 || chn==4 || chn==5))) && calib_en > 10000 && calib_en<8){
            bad_gamma=true;
        }
        
    }
    
    if( (ring_time-sector_time+200)>195 && (ring_time-sector_time+200)<215){
        ring_sector_time_gate=true;
    }
    
}


////////////////////////second loop fills histogram//////////////////
void second_loop(){
    
    Int_t entry_id;
    Int_t entry_size;
    Int_t count_rings_fired=0;
    Double_t time;
    Double_t en;
    Double_t calib_en;
    Double_t time_gate_lower=198;   //////////////////////////////////////////////////////////////////////remember to change normalisation for background subtraction
    Double_t time_gate_upper=211;
    Double_t time_gate_lower_La=198;
    Double_t time_gate_upper_La=211;
    Double_t time_gate_si_lower=198;   //used for efficiency
    Double_t time_gate_si_upper=211;   //used for efficiency
    Double_t diagonal_cut_width=2250;
    Double_t LaBr_uncalib_ds_gg=0;
    Double_t pg_gate_lower =203;
    Double_t pg_gate_upper =216;
    
    //Double_t time_gate_lower=0; //use when you don't want a time gate
    //Double_t time_gate_upper=1000;  //use when you don't want a time gate
    //ring_sector_time_gate=true; //set to true when you don't want a ring sector gate
    //ring_time=gamma_time;
    
    Double_t background_gate_lower=299;
    Double_t background_gate_upper=1101;
    Double_t fill_energy = 0.0;
    Double_t LaBr_uncalib_ds;
    Bool_t allready_filled;
    
    if (ring_sector_time_gate==true && abs(ring_en_for_hist - sector_en_for_hist)<diagonal_cut_width){
        ring_sector_time->Fill(sector_time-ring_time+200);
        //ring_en_sector_en->Fill(ring_en_for_hist, sector_en_for_hist);
    }
    
    
    for (entry_id = 0, entry_size = TRenergy.GetSize(); entry_id < entry_size; ++entry_id) {
        allready_filled = false;
        Int_t modchan=TRmodchan[entry_id];
        Int_t mod=modchan/16;
        Int_t chn=modchan%16;
        time=TRtime[entry_id];
        en=TRenergy[entry_id];
        if(mod==2&&(chn==0 || chn==1 || chn==2 || chn==3 || chn==4 || chn==5)){
            //for the LaBr we want to doppler shift first and then calibrate. For the Ge we want to calibrate, Add-back then doppler shift.
            //cout<< entry_id <<endl;
            calib_en=en;
        }
        else{
            calib_en=do_calibrate(en,mod,chn);
        }
        
        if( (mod==1 && (chn==9||chn==10||chn==11||chn==12||chn==13||chn==14||chn==15)) || (mod==2 && (chn==8||chn==9||chn==10||chn==11||chn==12||chn==13||chn==14||chn==15)) || (mod==3 && (chn==7||chn==8||chn==9||chn==10||chn==12||chn==13||chn==14||chn==15)) || mod==0 ){
            time_Si->Fill(time+200);
        }
        
        if( (mod==1 && (chn==8||chn==9||chn==10||chn==11||chn==12||chn==13||chn==14||chn==15))){
            rings_hit->Fill(chn-7);
        }
        if( (mod==2 && (chn==8||chn==9||chn==10||chn==11||chn==12||chn==13||chn==14||chn==15))){
            rings_hit->Fill(chn+1);
        }
        if( (mod==3 && (chn==7||chn==8||chn==9||chn==10))){
            rings_hit->Fill(chn+10);
        }
        if( (mod==3 && (chn==12||chn==13||chn==14||chn==15))){
            rings_hit->Fill(chn+9);
        }
        
        if( (mod==1&&(chn==0||chn==1||chn==2||chn==3))||(mod==3&&(chn==0||chn==1||chn==2||chn==3||chn==5)) || ( mod==2&&(chn==0||chn==1||chn==2||chn==3))){
            pg_time->Fill(time-ring_time+200);
        }
        
        if( (mod==1&&(chn==0||chn==1||chn==2||chn==3))||(mod==3&&(chn==0||chn==1||chn==2||chn==3)) ){
            // if ((time+200)>195 && (time+200)<210 && ring_sector_time_gate==true){
            time_Cl->Fill(time+200);
            //}
        }
        
        if(mod==3&&chn==5){
            time_Ge->Fill(time+200);
        }
        
        if( (mod==1&&(chn==0||chn==1||chn==2||chn==3))||(mod==3&&(chn==0||chn==1||chn==2||chn==3||chn==5)) ){
            //if ((time+200)>195 && (time+200)<210 && ring_sector_time_gate==true){
            time_Ge->Fill(time+200);
            //}
        }
        if( mod==2&&(chn==0||chn==1||chn==2||chn==3) ){
            time_La->Fill(time+200);
        }
        
        if (mod==3&&chn==13){
            one_ring_singles->Fill(calib_en);
            if((gamma_time+200)>time_gate_lower && (gamma_time+200)<time_gate_upper && abs(ring_en_for_hist - sector_en_for_hist)<diagonal_cut_width){
                one_ring_timegated->Fill(calib_en);
                if (abs(gamma_time-ring_time)<15){
                    one_ring_timegated_pggated->Fill(calib_en);
                }
            }
        }
        
        if( mod==2&&(chn==4||chn==5) ){
            //  if ((time+200)>195 && (time+200)<210 && ring_sector_time_gate==true){
            time_La_army->Fill(time+200);
            // }
        }
        
        hPixieADC[mod][chn]->Fill(calib_en); //This is not callibrated LaBr due to doppler being done before calibration
        
        if( (mod==1 && (chn==9||chn==10||chn==11||chn==12||chn==13||chn==14||chn==15)) || (mod==2 && (chn==8||chn==9||chn==10||chn==11||chn==12||chn==13||chn==14||chn==15)) || (mod==3 && (chn==7||chn==8||chn==9||chn==10||chn==12||chn==13||chn==14||chn==15))){
            //if((time-gamma_time+200)>199&&(time-gamma_time+200)<216&&time<25){
            rings_fired->Fill(ring_counter);
            //cout<<ring_counter<<endl;
            if (ring_counter==2){
                
                if (count_rings_fired==0){
                    ring_up->Fill(time);
                }
                if (count_rings_fired==1){
                    ring_middel->Fill(time);
                }
                if (count_rings_fired==2){
                    ring_down->Fill(time);
                }
                
                count_rings_fired++;
            }
            // }
        }
        if( mod==0 ){
            if((time-gamma_time+200)>199&&(time-gamma_time+200)<216){
                sectors_fired->Fill(sector_counter);
            }
        }
        
        ////////////////////////////energy sharing filling histogram/////////////////////////////////////
        if (energy_sharing_sector==true && chn==energy_sharing_sector_chn_use && mod==0&& abs(ring_en_for_hist - sector_en_for_hist)<diagonal_cut_width&&energy_sharing_sector_en>lower_lim&&energy_sharing_sector_en<upper_lim){
            
            if (abs(gamma_time-ring_time)<15 && (energy_share_use_time_sector+200)>time_gate_si_lower && (energy_share_use_time_sector+200)<time_gate_si_upper && ring_sector_time_gate==true&&calib_en>4000){
                hPixieADC_ab[mod][chn]->Fill(energy_sharing_sector_en);
                hPixieADC_pg[mod][chn]->Fill(energy_sharing_sector_en);
            }
            if (abs(gamma_time-ring_time)<15 && (energy_share_use_time_sector+200)>background_gate_lower && (energy_share_use_time_sector+200)<background_gate_upper && ring_sector_time_gate==true&&calib_en>4000){
                time_background[mod][chn]->Fill(energy_sharing_sector_en);
                time_background_not_dopler[mod][chn]->Fill(energy_sharing_sector_en);
            }
            allready_filled = true;
        }
        else if(energy_sharing_sector==true && chn==energy_sharing_sector_chn_ignore && mod==0){
            //do nothing
            allready_filled = true;
        }
        
        if (energy_sharing_ring==true && chn==energy_sharing_ring_chn_use && mod==energy_sharing_ring_mod_use && abs(ring_en_for_hist - sector_en_for_hist)<diagonal_cut_width&&energy_sharing_ring_en>lower_lim&&energy_sharing_ring_en<upper_lim){
            
            if (abs(gamma_time-ring_time)<15 && (energy_share_use_time_ring+200)>time_gate_si_lower && (energy_share_use_time_ring+200)<time_gate_si_upper && ring_sector_time_gate==true&&calib_en>4000){
                hPixieADC_pg[mod][chn]->Fill(energy_sharing_ring_en);
                hPixieADC_ab[mod][chn]->Fill(energy_sharing_ring_en);
            }
            if (abs(gamma_time-ring_time)<15 && (energy_share_use_time_ring+200)>background_gate_lower && (energy_share_use_time_ring+200)<background_gate_upper && ring_sector_time_gate==true&&calib_en>4000){
                time_background[mod][chn]->Fill(energy_sharing_ring_en);
                time_background_not_dopler[mod][chn]->Fill(energy_sharing_ring_en);
            }
            allready_filled = true;
        }
        else if(energy_sharing_ring==true && chn==energy_sharing_ring_chn_ignore && mod==energy_sharing_ring_mod_ignore){
            //do nothing
            //cout<<"spam2"<<endl;
            allready_filled = true;
        }
        /////////////////////////add_back histogram filling////////////////////////////////////////
        if (ab_event_cl1==true && chn==ab_use_crystal_cl1 && mod==1 && ring_sector_time_gate==true && abs(ring_en_for_hist - sector_en_for_hist)<diagonal_cut_width){
            beta=do_beta(beta_ring_en);
            fill_energy=do_doppler(ab_use_energy_cl1, mod, chn);
            if (abs(gamma_time-ring_time)<40 && (add_back_use_time1+200)>time_gate_lower && (add_back_use_time1+200)<time_gate_upper&&calib_en<10000){
                hPixieADC_pg[mod][chn]->Fill(fill_energy);
                hPixieADC_ab[mod][chn]->Fill(ab_use_energy_cl1);
                //ring_en_gamma_en->Fill(beta_ring_en, fill_energy);
                
            }
            if (abs(gamma_time-ring_time)<40 && (add_back_use_time1+200)>background_gate_lower && (add_back_use_time1+200)<background_gate_upper&&calib_en<10000){
                time_background[mod][chn]->Fill(fill_energy);
                time_background_not_dopler[mod][chn]->Fill(ab_use_energy_cl1);
                //ring_en_gamma_en_background->Fill(beta_ring_en, fill_energy);
            }
            
            allready_filled = true;
        }
        else if(ab_event_cl1==true && chn==ab_ignore_crystal_cl1 && mod==1){
            //do nothing
            allready_filled = true;
        }
        if (ab_event_cl3==true && chn==ab_use_crystal_cl3 && mod==3 && ring_sector_time_gate==true && abs(ring_en_for_hist - sector_en_for_hist)<diagonal_cut_width){
            beta=do_beta(beta_ring_en);
            fill_energy=do_doppler(ab_use_energy_cl3, mod, chn);
            
            if (abs(gamma_time-ring_time)<40 && (add_back_use_time3+200)>time_gate_lower && (add_back_use_time3+200)<time_gate_upper&&calib_en<10000 ){
                hPixieADC_pg[mod][chn]->Fill(fill_energy);
                hPixieADC_ab[mod][chn]->Fill(ab_use_energy_cl3);
                //ring_en_gamma_en->Fill(beta_ring_en, fill_energy);
            }
            if (abs(gamma_time-ring_time)<40 && (add_back_use_time3+200)>background_gate_lower && (add_back_use_time3+200)<background_gate_upper&&calib_en<10000){
                time_background[mod][chn]->Fill(fill_energy);
                time_background_not_dopler[mod][chn]->Fill(ab_use_energy_cl3);
                //ring_en_gamma_en_background->Fill(beta_ring_en, fill_energy);
            }
            
            
            allready_filled = true;
        }
        else if(ab_event_cl3==true && chn==ab_ignore_crystal_cl3 && mod==3){
            //do nothing
            allready_filled = true;
        }
        
        /////////////////////////everything else histogram filling////////////////////////////////////////
        
        if(allready_filled==false && calib_en>20 && ring_sector_time_gate==true && abs(ring_en_for_hist - sector_en_for_hist)<diagonal_cut_width){
            
            if (energy_sharing_ring==true){
                ring_mod = energy_sharing_ring_mod_use;
                ring_chn = energy_sharing_ring_chn_use;
            }
            if (energy_sharing_sector==true){
                sector_chn=energy_sharing_sector_chn_use;
                sector_mod=0;
            }
            
            if(mod == energy_sharing_ring_mod_ignore && chn == energy_sharing_ring_chn_ignore){
                //ignore
            }
            if ((time-ring_time+200)>203&&(time-ring_time+200)<216){
                
                if(mod==2&&(chn==0 || chn==1 || chn==2 || chn==3)&&calib_en<10000 && (time+200)> time_gate_lower_La && (time+200)<time_gate_upper_La){
                    //for the LaBr we want to doppler shift first and then calibrate. For the Ge we want to calibrate, Add-back then doppler shift.
                    beta=do_beta(beta_ring_en);
                    LaBr_uncalib_ds=do_doppler(calib_en, mod, chn);
                    fill_energy=do_calibrate(LaBr_uncalib_ds,mod,chn);
                    calib_en=do_calibrate(en,mod,chn);
                    ring_en_gamma_en->Fill(beta_ring_en, fill_energy);
                    hPixieADC_pg[mod][chn]->Fill(fill_energy);
                    hPixieADC_ab[mod][chn]->Fill(calib_en);
                    //ring_en_sector_en->Fill(ring_en_for_hist, sector_en_for_hist);
                }
                if(((mod==1&&(chn==0 || chn==1 || chn==2 || chn==3))||(mod==3&&(chn==0 || chn==1 || chn==2 || chn==3 || chn==5)))&& calib_en<10000   && (time+200)>time_gate_lower && (time+200)<time_gate_upper){
                    beta=do_beta(beta_ring_en);
                    fill_energy=do_doppler(calib_en, mod, chn);
                    //ring_en_gamma_en->Fill(beta_ring_en, fill_energy);
                    hPixieADC_pg[mod][chn]->Fill(fill_energy);
                    hPixieADC_ab[mod][chn]->Fill(calib_en);
                    
                    //ring_en_sector_en->Fill(ring_en_for_hist, sector_en_for_hist);
                }
                else if((time+200)>time_gate_si_lower && (time+200)<time_gate_si_upper && calib_en>10){ //i.e. if it is the silicon
                    //fill_energy=calib_en;
                    hPixieADC_pg[mod][chn]->Fill(calib_en);
                    hPixieADC_ab[mod][chn]->Fill(calib_en);
                    ring_en_sector_en->Fill(ring_en_for_hist, sector_en_for_hist);
                    /*if (sector_en_for_hist<lower_lim||sector_en_for_hist>upper_lim||ring_en_for_hist<lower_lim||ring_en_for_hist>upper_lim){
                     cout<<"sector: "<<sector_en_for_hist<<" ring: "<<ring_en_for_hist<<endl;
                     }*/
                }
                
                
                //time background
                if(mod==2&&(chn==0 || chn==1 || chn==2 || chn==3)&&calib_en<10000 && (time+200)>background_gate_lower && (time+200)<background_gate_upper){
                    //for the LaBr we want to doppler shift first and then calibrate. For the Ge we want to calibrate, Add-back then doppler shift.
                    beta=do_beta(beta_ring_en);
                    LaBr_uncalib_ds=do_doppler(calib_en, mod, chn);
                    fill_energy=do_calibrate(LaBr_uncalib_ds,mod,chn);
                    ring_en_gamma_en_background->Fill(beta_ring_en, fill_energy);
                    calib_en=do_calibrate(calib_en,mod,chn);
                    time_background[mod][chn]->Fill(fill_energy);
                    time_background_not_dopler[mod][chn]->Fill(calib_en);
                }
                else if(((mod==1&&(chn==0 || chn==1 || chn==2 || chn==3))||(mod==3&&(chn==0 || chn==1 || chn==2 || chn==3 || chn==5)))&&calib_en<10000 && (time+200)>background_gate_lower && (time+200)<background_gate_upper){
                    beta=do_beta(beta_ring_en);
                    fill_energy=do_doppler(calib_en, mod, chn);
                    //ring_en_gamma_en_background->Fill(beta_ring_en, fill_energy);
                    time_background[mod][chn]->Fill(fill_energy);
                    time_background_not_dopler[mod][chn]->Fill(calib_en);
                    
                }
                else if((time+200)>background_gate_lower && (time+200)<background_gate_upper&&calib_en>4000){ //i.e. if it is the silicon
                    //fill_energy=calib_en;
                    time_background[mod][chn]->Fill(calib_en);
                    time_background_not_dopler[mod][chn]->Fill(calib_en);
                }
                
            }
            
            else if((time-ring_time+200)>(pg_gate_upper+50)&&(time-ring_time+200)<(pg_gate_upper+900)){ //gamma-particle background
                
                if(mod==2&&(chn==0 || chn==1 || chn==2 || chn==3)&&calib_en<10000 && (time+200)> time_gate_lower_La && (time+200)<time_gate_upper_La){
                    //for the LaBr we want to doppler shift first and then calibrate. For the Ge we want to calibrate, Add-back then doppler shift.
                    beta=do_beta(beta_ring_en);
                    LaBr_uncalib_ds=do_doppler(calib_en, mod, chn);
                    fill_energy=do_calibrate(LaBr_uncalib_ds,mod,chn);
                    ring_en_gamma_en_background->Fill(beta_ring_en, fill_energy);
                    calib_en=do_calibrate(calib_en,mod,chn);
                    time_background_pg[mod][chn]->Fill(fill_energy);
                    time_background_not_dopler_pg[mod][chn]->Fill(calib_en);
                    
                    
                    //ring_en_sector_en->Fill(ring_en_for_hist, sector_en_for_hist);
                }
                else if(((mod==1&&(chn==0 || chn==1 || chn==2 || chn==3))||(mod==3&&(chn==0 || chn==1 || chn==2 || chn==3 || chn==5)))&& calib_en<10000 && (time+200)>time_gate_lower && (time+200)<time_gate_upper){
                    beta=do_beta(beta_ring_en);
                    fill_energy=do_doppler(calib_en, mod, chn);
                    //ring_en_gamma_en_background->Fill(beta_ring_en, fill_energy);
                    time_background_pg[mod][chn]->Fill(fill_energy);
                    time_background_not_dopler_pg[mod][chn]->Fill(calib_en);
                    //ring_en_sector_en->Fill(ring_en_for_hist, sector_en_for_hist);
                }
                else if((time+200)>time_gate_si_lower && (time+200)<time_gate_si_upper && calib_en>10){ //i.e. if it is the silicon
                    //fill_energy=calib_en;
                    time_background_pg[mod][chn]->Fill(calib_en);
                    time_background_not_dopler_pg[mod][chn]->Fill(calib_en);
                }
                
                
                //time background
                if(mod==2&&(chn==0 || chn==1 || chn==2 || chn==3)&&calib_en<10000 && (time+200)>background_gate_lower && (time+200)<background_gate_upper){
                    //for the LaBr we want to doppler shift first and then calibrate. For the Ge we want to calibrate, Add-back then doppler shift.
                    beta=do_beta(beta_ring_en);
                    LaBr_uncalib_ds=do_doppler(calib_en, mod, chn);
                    fill_energy=do_calibrate(LaBr_uncalib_ds,mod,chn);
                    ring_en_gamma_en_background->Fill(beta_ring_en, fill_energy);
                    calib_en=do_calibrate(calib_en,mod,chn);
                    time_background_pg[mod][chn]->Fill(fill_energy);
                    time_background_not_dopler_pg[mod][chn]->Fill(calib_en);
                }
                else if(((mod==1&&(chn==0 || chn==1 || chn==2 || chn==3))||(mod==3&&(chn==0 || chn==1 || chn==2 || chn==3 || chn==5)))&&calib_en<10000 && (time+200)>background_gate_lower && (time+200)<background_gate_upper){
                    beta=do_beta(beta_ring_en);
                    fill_energy=do_doppler(calib_en, mod, chn);
                    //ring_en_gamma_en_background->Fill(beta_ring_en, fill_energy);
                    time_background_pg[mod][chn]->Fill(fill_energy);
                    time_background_not_dopler_pg[mod][chn]->Fill(calib_en);
                    
                }
                else if((time+200)>background_gate_lower && (time+200)<background_gate_upper&&calib_en>4000){ //i.e. if it is the silicon
                    //fill_energy=calib_en;
                    time_background_pg[mod][chn]->Fill(calib_en);
                    time_background_not_dopler_pg[mod][chn]->Fill(calib_en);
                }
                
            }
        }
        
    }
}

void read_calibrate_file() {
    Double_t a = 0.0;
    Double_t b = 1.0;
    Double_t c = 1.0;
    Double_t d = 0.0;
    Double_t e = 0.0;
    Double_t f = 0.0;
    Double_t g = 0.0;
    FILE *calib = fopen("Calibration_before_exp_14C.txt","r");
    //FILE *calib = fopen("Calibration_empty.txt","r");
    for(int m=0;m<4;m++){
        for(int n=0;n<16;n++){ //The format is: offset(a) slope^1(b) slope(c)^2 phi(d) theta(e)
            fscanf(calib, "%lf %lf %lf %lf %lf %lf %lf", &a, &b, &c, &f, &g, &d, &e);
            calib_a[m][n]=a;
            calib_b[m][n]=b;
            calib_c[m][n]=c;
            calib_d[m][n]=f;
            calib_e[m][n]=g;
            phi[m][n]=d;
            theta[m][n]=e;
        }
    }
    
}

Double_t do_calibrate(Double_t raw_energy, Int_t module, Int_t channel){
    Double_t random_num;
    Double_t cal_energy;
    Double_t offset = calib_a[module][channel];
    Double_t slope1 = calib_b[module][channel];
    Double_t slope2 = calib_c[module][channel];
    Double_t slope3 = calib_d[module][channel];
    Double_t slope4 = calib_e[module][channel];
    random_num=r->Uniform(-0.5,0.5);
    cal_energy=(raw_energy+random_num)*(raw_energy+random_num)*(raw_energy+random_num)*(raw_energy+random_num)*slope4+(raw_energy+random_num)*(raw_energy+random_num)*(raw_energy+random_num)*slope3+(raw_energy+random_num)*(raw_energy+random_num)*slope2 + (raw_energy+random_num)*slope1 + offset;
    //cal_energy=raw_energy+random_num;
    //if (module==2&&channel==14){
    //    cout << raw_energy << " " << cal_energy << endl;
    //}
    return cal_energy;
}


Double_t do_beta(Double_t ring_en){
    Double_t beta;
    beta=sqrt(((ring_en/1000.0)*2.0)/(931.5*14.0));
    //beta=0.071;
    /*if (beta==0){
     cout<<ring_en<<endl;
     }*/
    //cout<<ring_en<< " " << beta<<endl;
    return beta;
}


Double_t do_doppler(Double_t en, Int_t gamma_mod, Int_t gamma_chn){
    Double_t dop_cor_en;
    Double_t cos_alpha;
    cos_alpha=sin(theta[ring_mod][ring_chn])*sin(theta[gamma_mod][gamma_chn])*cos(phi[sector_mod][sector_chn]-phi[gamma_mod][gamma_chn]) + cos(theta[ring_mod][ring_chn])*cos(theta[gamma_mod][gamma_chn]);
    //cos_alpha=cos(theta[gamma_mod][gamma_chn]-theta[ring_mod][ring_chn]);
    dop_cor_en=en*(1.0+beta*cos_alpha)/(sqrt(1-beta*beta));
    /*if (gamma_mod==2&&gamma_chn==2&&en>6000){
     cout << "cos alpha " << cos_alpha << " energy before " << en << " doppler shifted energy " << dop_cor_en << " theta_p " << theta[ring_mod][ring_chn]  << " phi_p " << phi[sector_mod][sector_chn] << " phi_gamma " << phi[gamma_mod][gamma_chn]<< " theta_gamma " << theta[gamma_mod][gamma_chn] << " beta " << beta << endl;
     
     }*/
    /*if (en>5800 && print_army_count<50 && gamma_mod == 2 && (gamma_chn==0) ){
     cout << "Energy before: " << en << " DS energy: " << dop_cor_en << " Alpha: " << cos_alpha <<  " Theta_p " << theta[ring_mod][ring_chn]  << " Phi_p " << phi[sector_mod][sector_chn] << " theta_gamma " << theta[gamma_mod][gamma_chn]<< " phi_gamma " << phi[gamma_mod][gamma_chn] << " Beta " << beta << " Det mod " << gamma_mod << " chn, " << gamma_chn << endl;
     print_army_count++;
     }
     if (en>5800 && print_oslo_count<50 && gamma_mod == 2 && (gamma_chn==1) ){
     cout << "Energy before: " << en << " DS energy: " << dop_cor_en << " Alpha: " << cos_alpha <<  " Theta_p " << theta[ring_mod][ring_chn]  << " Phi_p " << phi[sector_mod][sector_chn] << " theta_gamma " << theta[gamma_mod][gamma_chn]<< " phi_gamma " << phi[gamma_mod][gamma_chn] << " Beta " << beta << " Det mod " << gamma_mod << " chn, " << gamma_chn << endl;
     print_oslo_count++;
     }
     if (en>5800 && print_clover_count<50 && (gamma_mod == 2 && gamma_chn == 5)){
     cout << "Energy before: " << en << " DS energy: " << dop_cor_en << " Alpha: " << cos_alpha <<  " Theta_p " << theta[ring_mod][ring_chn]  << " Phi_p " << phi[sector_mod][sector_chn] << " theta_gamma " << theta[gamma_mod][gamma_chn]<< " phi_gamma " << phi[gamma_mod][gamma_chn] << " Beta " << beta << " Det mod " << gamma_mod << " chn, " << gamma_chn << endl;
     print_clover_count++;
     }*/
    return dop_cor_en;
}


void init_histos() {
    //calibrated
    /*Int_t gamma_bins=10000;
     Int_t gamma_min=0;
     Int_t gamma_max=9999;
     Int_t particle_bins=36000;
     Int_t particle_min=0;
     Int_t particle_max=35999;*/
    
    //uncalibrated
    Int_t gamma_bins=32769;
    Int_t gamma_min=0;
    Int_t gamma_max=32768;
    Int_t particle_bins=48000;
    Int_t particle_min=0;
    Int_t particle_max=47999;
    
    ring_up= new TH1F("ring_up","",particle_bins,particle_min,particle_max);
    ring_middel= new TH1F("ring_middel","",particle_bins,particle_min,particle_max);
    ring_down= new TH1F("ring_down","",particle_bins,particle_min,particle_max);
    time_Ge = new TH1F("time_Ge_single","",1501,0,1500);
    time_Cl = new TH1F("time_Cl","",1501,0,1500);
    time_La = new TH1F("time_La","",1501,0,1500);
    time_La_army = new TH1F("time_La_army","",1501,0,1500);
    ring_sector_time = new TH1F("time_ring_sector","",1501,0,1500);
    pg_time = new TH1F("time_particle_gamma","",1501,0,1500);
    time_Si = new TH1F("time_Si","",1501,0,1500);
    Clover1 = new TH1F("Clover_1","",particle_bins,particle_min,particle_max);
    Clover3 = new TH1F("Clover_3","",particle_bins,particle_min,particle_max);
    Clover1ab = new TH1F("Clover_1ab","",particle_bins,particle_min,particle_max);
    Clover3ab = new TH1F("Clover_3ab","",particle_bins,particle_min,particle_max);
    one_ring_singles = new TH1F("one_ring_singles","",particle_bins,particle_min,particle_max);
    one_ring_timegated = new TH1F("one_ring_timegated","",particle_bins,particle_min,particle_max);
    one_ring_timegated_pggated = new TH1F("one_ring_timegated_pggated","",particle_bins,particle_min,particle_max);
    rings_hit = new TH1F("rings_hit","",26,0,25);
    ring_en_sector_en = new TH2F("ring_en_sector_en","",1000,particle_min,48000,1000,particle_min,48000);
    ring_en_gamma_en = new TH2F("ring_en_gamma_en","",1000,particle_min,particle_max,1000,particle_min,9000);
    ring_en_gamma_en_background = new TH2F("ring_en_gamma_en_background","",1000,particle_min,particle_max,1000,particle_min,9000);
    all_detectors_complete= new TH1F("all_detectors_complete","",particle_bins,particle_min,particle_max);
    all_detectors_ab= new TH1F("all_detectors_ab","",particle_bins,particle_min,particle_max);
    all_Ge_ab= new TH1F("all_Ge_ab","",particle_bins,particle_min,particle_max);
    all_LaBr_ab= new TH1F("all_LaBr_ab","",particle_bins,particle_min,particle_max);
    all_Ge= new TH1F("all_Ge","",particle_bins,particle_min,particle_max);
    all_LaBr= new TH1F("all_LaBr","",particle_bins,particle_min,particle_max);
    rings_fired = new TH1I("rings_fired","",10,0,10);
    sectors_fired = new TH1I("sectors_fired","",10,0,10);
    
    for(Int_t j=0;j<4;j++){
        
        for(Int_t i=0;i<16;i++){
            hPixieADC_pg[j][i] = new TH1F(Form("ADC_complete_%d_ch%02d",j,i),"",particle_bins,particle_min,particle_max);
        }
    }
    for(Int_t j=0;j<4;j++){
        
        for(Int_t i=0;i<16;i++){
            hPixieADC_ab[j][i] = new TH1F(Form("ADC_ab%d_ch%02d",j,i),"",particle_bins,particle_min,particle_max);
        }
    }
    for(Int_t j=0;j<4;j++){
        
        for(Int_t i=0;i<16;i++){
            hPixieADC[j][i] = new TH1F(Form("ADC%d_ch%02d",j,i),"",particle_bins,particle_min,particle_max);
        }
    }
    for(Int_t j=0;j<4;j++){
        
        for(Int_t i=0;i<16;i++){
            time_background[j][i] = new TH1F(Form("time_background%d_ch%02d",j,i),"",particle_bins,particle_min,particle_max);
        }
    }
    
    for(Int_t j=0;j<4;j++){
        
        for(Int_t i=0;i<16;i++){
            time_background_not_dopler[j][i] = new TH1F(Form("time_background_not_dopler%d_ch%02d",j,i),"",particle_bins,particle_min,particle_max);
        }
    }
    for(Int_t j=0;j<4;j++){
        
        for(Int_t i=0;i<16;i++){
            time_background_pg[j][i] = new TH1F(Form("time_background%d_ch%02d_pg",j,i),"",particle_bins,particle_min,particle_max);
        }
    }
    
    for(Int_t j=0;j<4;j++){
        
        for(Int_t i=0;i<16;i++){
            time_background_not_dopler_pg[j][i] = new TH1F(Form("time_background_not_dopler%d_ch%02d_pg",j,i),"",particle_bins,particle_min,particle_max);
        }
    }
    /*for(Int_t j=0;j<4;j++){
     
     for(Int_t i=0;i<16;i++){
     hPixieADC_eff[j][i] = new TH1F(Form("ADC_eff_%d_ch%02d",j,i),"",particle_bins,particle_min,particle_max);
     }
     }
     
     for(Int_t j=0;j<4;j++){
     
     for(Int_t i=0;i<16;i++){
     time_background_eff[j][i] = new TH1F(Form("time_background_eff%d_ch%02d",j,i),"",particle_bins,particle_min,particle_max);
     }
     }*/
    
    
}

void make_chain(){
    
    //////14C/////////////
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0675*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0676*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0677*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0678*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0679*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-068*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-069*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-070*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-071*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-072*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-073*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-074*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-075*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-076*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-077*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-078*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0790*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0791*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0809*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0810*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0811*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0812*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0813*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0814*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0815*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0816*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0817*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0818*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0822*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0823*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0824*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0825*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0826*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0827*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0828*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0829*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-083*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-084*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-085*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-086*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-087*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-088*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-089*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0901*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0902*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0903*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0904*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0905*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0906*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0907*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0908*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0910*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0911*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0912*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0913*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0914*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0915*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0916*.root");
    
    //DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0689-02*.root");   //one for testing (900mb)
    
}
